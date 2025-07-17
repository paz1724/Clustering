classdef cFCM < handle
    %FCMCLUSTERER  Fuzzy C-Means clustering with auto-selection and outlier removal
    %
    % Tunable Properties:
    %   numClusters           - # clusters (c). If empty & autoSelect=true, determined via optimize.
    %   autoSelect            - logical; if true, pick numClusters automatically.
    %   clusterRange          - vector of candidate c values when autoSelect=true.
    %   minClusterSize        - enforce at least this many points per cluster
    %   validityIndex         - 'XB' (Xie–Beni) or 'PC' (Partition Coeff).
    %   fuzzifier             - fuzziness exponent m (>1).
    %   maxIter               - max iterations for FCM.
    %   tol                   - center‐shift tolerance.
    %   distanceMetric        - @(X,V) → N×c distances between data X and centers V.
    %   outlierMethod         - 'iqr' or 'percentile'.
    %   outlierFraction       - fraction to keep if 'percentile'.
    %   outlierIQRFactor      - IQR multiplier if 'iqr'.
    %   iqrLowerPercentile    - lower percentile for IQR (default 25).
    %   iqrUpperPercentile    - upper percentile for IQR (default 75).
    %
    % After fit(X), these are set:
    %   centers     - c×d cluster centers.
    %   labels      - N×1 hard labels (1…c).
    %   cleanLabels - N×1 labels after outlier pruning (0 = removed).
    %   outlierIdx  - N×1 logical mask of removed points.

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~ Properties ~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    properties
        numClusters        = []
        autoSelect         = false
        clusterRange       = 1:6
        minClusterSize     = 3;
        validityIndex      = 'XB'
        fuzzifier          = 2.0
        maxIter            = 100
        tol                = 1e-5
        distanceMetric     = @(X,V) pdist2(X,V)
        outlierMethod      = 'iqr'
        outlierFraction    = 0.90
        outlierIQRFactor   = 1.5
        iqrLowerPercentile = 25
        iqrUpperPercentile = 75
        shapePenaltyWeight = 0          % =0 turns off shape-penalty
        costTh             = 0.5
        shapePenaltyMethod = 'axisRatio';% 'axisRatio' or 'volume'
    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~ Private Properties ~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    properties (SetAccess=private)
        centers
        labels
        cleanLabels
        outlierIdx
    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~ Methods ~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods
        %% ~~~~~~~~~~~~~~~ FCMClusterer ~~~~~~~~~~~~~~~ %%
        function obj = cFCM(varargin)
            % Allow name/value assignments
            if mod(numel(varargin),2)~=0
                error('Name/value pairs required.');
            end
            for k=1:2:numel(varargin)
                prop = varargin{k}; val = varargin{k+1};
                if isprop(obj,prop)
                    obj.(prop) = val;
                else
                    error('Unknown property "%s".', prop);
                end
            end
        end

        %% ~~~~~~~~~~~~~~~ fit ~~~~~~~~~~~~~~~ %%
        function Fit(obj, X)

            %% Produce Num Clusters
            if obj.autoSelect
                obj.numClusters = obj.Select_optimal_C(X);
                fprintf('Auto‐selected numClusters = %d\n', obj.numClusters);
            elseif isempty(obj.numClusters)
                error('numClusters empty & autoSelect=false.');
            end

            %% Run
            [obj.labels, obj.centers] = obj.Run(obj.numClusters, X);

            %% Prune_outliers
            obj.Prune_outliers(X);

            %% Remove_outlier_classes
            obj.Remove_outlier_classes();

        end

        %% ~~~~~~~~~~~~~~~ plotResults ~~~~~~~~~~~~~~~ %%
        function Plot_results(obj, X, gtMeans, seed)
            %PLOTRESULTS  Show raw vs. pruned clustering
            if isempty(obj.labels)
                error('Call fit(X) first.');
            end
            c = obj.numClusters;
            cols = lines(c);

            figure;

            %% Plot Raw
            subplot(1,2,1);
            hold on;
            grid minor;

            %% Scatter GT
            scatter(gtMeans(:,1), gtMeans(:,2), 100, 'k', 'Marker', 's','LineWidth',1.5);

            for k=1:c
                scatter(X(obj.labels==k,1), X(obj.labels==k,2), 36, cols(k,:), 'filled');
                cnt = sum(obj.labels==k);
                if cnt==0
                    marker = 'x';
                else
                    marker = '*';
                end
                scatter(obj.centers(k,1), obj.centers(k,2), 100, cols(k,:), 'Marker', marker,'LineWidth',1.5);
            end
            legend([{'GT'},arrayfun(@(k)sprintf('C%d',ceil(k/2)),1:2*c,'uni',false),{'Centers'}], ...
                'Location','BestOutside');
            hold off;

            % build multiline title: heading + centers
            info = cell(c+1,1);
            info{1} = 'FCM (raw)';
            for k=1:c
                cnt = sum(obj.labels==k);
                ctr = obj.centers(k,:);
                info{k+1} = sprintf('C%d Center = (%.2f, %.2f) --> n = %d', k, ctr(1), ctr(2), cnt);
            end
            title(info, 'Interpreter','none');

            %% Plot Pruned
            subplot(1,2,2);
            hold on;
            grid minor;

            %% Scatter GT
            scatter(gtMeans(:,1), gtMeans(:,2), 100, 'k', 'Marker', 's','LineWidth',1.5);
 
            for k=1:c
                scatter(X(obj.cleanLabels==k,1), X(obj.cleanLabels==k,2), 36, cols(k,:), 'filled');

                cnt = sum(obj.cleanLabels==k);
                if cnt==0
                    marker = 'x';
                else
                    marker = '*';
                end
                scatter(obj.centers(k,1), obj.centers(k,2), 100, cols(k,:), 'Marker', marker,'LineWidth',1.5);
            end
            scatter(X(obj.outlierIdx,1), X(obj.outlierIdx,2), 36, [.5 .5 .5], 'x');

            legend([{'GT'},arrayfun(@(k)sprintf('C%d',ceil(k/2)),1:2*c,'uni',false),{'Outliers','Centers'}], ...
                'Location','BestOutside');
            hold off;

            % Cluster info: center and count per cluster, one row each
            info = cell(1,1);
            for k=1:c
                cnt = sum(obj.cleanLabels==k);
                ctr = obj.centers(k,:);
                info{k,1} = sprintf('C%d: Center = (%.2f,%.2f) --> n = %d', k, ctr(1), ctr(2), cnt);
            end
            titleList = vertcat({'FCM (Pruned)';"Seed: "+seed},info);
            title(titleList, 'Interpreter','none');

            plotbrowser('on');
        end

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~ Private Methods ~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods (Access=private)
        %% ~~~~~~~~~~~~~~~ Run ~~~~~~~~~~~~~~~ %%
        function [labels, centers] = Run(obj, c, X)
            % Vectorized FCM: hard labels + centers
            [N,d] = size(X);
            p     = 2/(obj.fuzzifier-1);
            U     = rand(c,N);
            U     = bsxfun(@rdivide, U, sum(U,1));
            v     = (U.^obj.fuzzifier * X) ./ (sum(U.^obj.fuzzifier,2)*ones(1,d));
            vOld  = v;
            delta = inf;
            it    = 0;
            while it<obj.maxIter && delta>obj.tol
                Dmat  = obj.distanceMetric(X, v);      % N×c
                DP    = Dmat.^(-p);                    % N×c
                U     = bsxfun(@rdivide, DP', sum(DP,2)'); % c×N
                vnew  = (U.^obj.fuzzifier * X) ./ (sum(U.^obj.fuzzifier,2)*ones(1,d));
                delta = max(abs(vnew(:)-vOld(:)));
                vOld  = vnew;
                v     = vnew;
                it    = it+1;
            end
            [powers, labs] = max(U,[],1);

            labs(powers<obj.costTh) = 0;
            labels  = labs(:);
            centers = v;
        end

        %% ~~~~~~~~~~~~~~~ Select_optimal_C ~~~~~~~~~~~~~~~ %%
        function bestC = Select_optimal_C(obj, X)
            % Brute‐force search over clusterRange
            bestVal = Inf;
            bestC   = obj.clusterRange(1);
            for c   = obj.clusterRange
                val = obj.Validity_cost(c, X);
                if val < bestVal
                    bestVal = val;
                    bestC   = c;
                end
            end
        end

        %% ~~~~~~~~~~~~~~~ validityCost ~~~~~~~~~~~~~~~ %%
        function val = Validity_cost(obj, c, X)
            % Compute validity index for integer c (including c=1)
            [U,V] = obj.Compute_fuzzy(c, X);
            D     = obj.distanceMetric(X, V)';  % c×N
            Um    = U.^obj.fuzzifier;
            N     = size(X,1);
            switch upper(obj.validityIndex)
                case 'XB'
                    num = sum(sum( Um .* (D.^2) ));
                    if c == 1
                        % for a single cluster, define XB as average within‐cluster SSE
                        val = num / N;
                    else
                        d2                  = pdist2(V,V).^2;
                        d2(logical(eye(c))) = Inf;
                        denom               = N * min(d2(:));
                        val                 = num / denom;
                    end
                case 'PC'
                    val = sum(sum(U.^2)) / N;
                otherwise
                    error('Unknown validityIndex "%s".', obj.validityIndex);
            end

            if obj.shapePenaltyWeight>0
                sp = obj.Compute_Shape_Penalty(X, V, U);  % c×1 penalties
                val = val + obj.shapePenaltyWeight * mean(sp);
            end
        end

        %% ~~~~~~~~~~~~~~~ Compute_fuzzy ~~~~~~~~~~~~~~~ %%
        function [U, V] = Compute_fuzzy(obj, c, X)
            % Return full membership U (c×N) and centers V (c×d)
            [N,d] = size(X);
            p     = 2/(obj.fuzzifier-1);
            U     = rand(c,N);
            U     = bsxfun(@rdivide, U, sum(U,1));
            delta = inf;
            it    = 0;
            while it<obj.maxIter && delta>obj.tol
                Um    = U.^obj.fuzzifier;
                V     = (Um * X) ./ (sum(Um,2)*ones(1,d));
                Dmat  = obj.distanceMetric(X, V);
                DP    = Dmat.^(-p);
                Unew  = bsxfun(@rdivide, DP', sum(DP,2)');
                delta = max(abs(Unew(:)-U(:)));
                U     = Unew;
                it    = it+1;
            end
        end

        %% ~~~~~~~~~~~~~~~ Prune_outliers ~~~~~~~~~~~~~~~ %%
        function Prune_outliers(obj, X)
            %PRUNEOUTLIERS  Remove per-cluster outliers by chosen rule
            if isempty(obj.labels)
                error('Call fit(X) first.');
            end
            N      = size(X,1); K = obj.numClusters;
            clean  = obj.labels;
            outidx = false(N,1);
            for k=1:K
                idx = find(obj.labels==k);
                if isempty(idx), continue; end
                d = sqrt(sum((X(idx,:)-obj.centers(k,:)).^2,2));
                switch lower(obj.outlierMethod)
                    case 'iqr'
                        Q1  = prctile(d, obj.iqrLowerPercentile);
                        Q3  = prctile(d, obj.iqrUpperPercentile);
                        thr = Q3 + obj.outlierIQRFactor*(Q3 - Q1);
                    case 'percentile'
                        thr = prctile(d, obj.outlierFraction*100);
                    otherwise
                        error('Unknown outlierMethod "%s".', obj.outlierMethod);
                end
                bad = d > thr;
                outidx(idx(bad)) = true;
                clean(idx(bad))  = 0;
            end
            obj.cleanLabels = clean;
            obj.outlierIdx  = outidx;
        end

        %% ~~~~~~~~~~~~~~~ Remove_outlier_classes ~~~~~~~~~~~~~~~ %%
        function [] = Remove_outlier_classes(obj)
            counts = histcounts(obj.labels, 0.5:1:(obj.numClusters+0.5));
            badClusters = find(counts < obj.minClusterSize);
            if ~isempty(badClusters)
                % any point whose label is in badClusters becomes noise:
                mask = ismember(obj.labels, badClusters);
                % obj.labels(mask)      = 0;
                obj.cleanLabels(mask,1) = 0;
                obj.outlierIdx(mask,1)  = true;
            end
        end

        %% ~~~~~~~~~~~~~~~ Compute_Shape_Penalty ~~~~~~~~~~~~~~~ %%
        function p = Compute_Shape_Penalty(obj, X, V, U)
            % compute penalty for each cluster i
            % X : N×d data
            % V : c×d centers
            % U : c×N memberships (for hard labels just do [~,L]=max(U); idx=find(L==i))
            c = size(V,1);
            p = zeros(c,1);
            % if you prefer hard labels:
            [~, L] = max(U,[],1);
            for i=1:c
                idx = find(L==i);
                if numel(idx)<2
                    p(i) = 0;
                    continue;
                end
                C = cov(X(idx,:));
                switch lower(obj.shapePenaltyMethod)
                    case 'axisratio'
                        e = eig(C);
                        p(i) = sqrt(max(e)/(min(e)+1e-5));
                    case 'volume'
                        p(i) = det(C);
                    otherwise
                        error('Unknown shapePenaltyMethod "%s".',obj.shapePenaltyMethod);
                end
            end
        end

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~ Static Methods ~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods (Static)
        %% ~~~~~~~~~~~~~~~ Calc_GDOP ~~~~~~~~~~~~~~~ %%
        function gdops = Calc_GDOP(centers, anchors)
            %MYGDOPCALCULATOR  Compute 2D GDOP for each center given anchor positions.
            %   centers : c×2
            %   anchors : m×2
            %   gdops    : c×1

            [c,~]      = size(centers);
            [m,~]      = size(anchors);
            gdops      = zeros(c,1);

            for i = 1:c
                % Vector from center i to each anchor
                dx     = anchors(:,1) - centers(i,1);
                dy     = anchors(:,2) - centers(i,2);
                r      = sqrt(dx.^2 + dy.^2);     % ranges
                % build geometry matrix H (m×2) of direction cosines
                H      = [dx./r, dy./r];
                % covariance of a least-squares fix is (H′H)⁻¹
                Q      = inv(H'*H);
                % GDOP is sqrt(trace(Q))
                gdops(i) = sqrt( trace(Q) );
            end

        end
    end
end
