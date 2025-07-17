classdef (Abstract) cClusterer < handle
    %cCLUSTERER  Abstract base for clustering algorithms
    %
    % Shared Properties (populated by Fit):
    %   centers     – c×d cluster center coordinates
    %   labels      – N×1 integer labels (0 = noise/unassigned)
    %   cleanLabels – N×1 labels after outlier pruning (0 = removed)
    %   outlierIdx  – N×1 logical mask of outlier points
    %
    % Subclasses must implement:
    %   Fit(obj, X)
    %       – perform clustering on data X (N×d)
    %       – populate the four properties above
    %
    % Provided Method:
    %   Plot_results(obj, X, sPlot)
    %       – visualize raw vs. pruned clustering
    %       – sPlot is a struct with fields:
    %           .gtMeans – g×d ground‐truth centers (optional)
    %           .seed    – scalar seed for title
    %           .algMode – string algorithm name for title

    properties
        centers
        labels
        cleanLabels
        outlierIdx
    end

    methods (Abstract)
        Fit(obj, X)
    end

    methods
        %% ~~~~~~~~~~~~~~~ Plot_results ~~~~~~~~~~~~~~~ %%
        function Plot_results(obj, X, sPlot)
            % PLOTRESULTS  Show raw vs. pruned clustering
            if isempty(obj.labels)
                error('Call Fit(X) first.');
            end
            c       = size(obj.centers,1);
            cols    = lines(c);
            gtMeans = sPlot.gtMeans;
            seed    = sPlot.seed;
            algMode = sPlot.algMode;

            figure;

            %% Plot Raw
            subplot(1,2,1);
            hold on;
            grid minor;

            %% Scatter GT
            scatter(gtMeans(:,1), gtMeans(:,2), 100, 'k', 's','LineWidth',1.5);

            for k=1:c
                scatter(X(obj.labels==k,1), X(obj.labels==k,2), 36, cols(k,:), 'filled');
                cnt = sum(obj.labels==k);
                marker = ternary(cnt>0, '*', 'x');
                scatter(obj.centers(k,1), obj.centers(k,2), 100, cols(k,:), ...
                        'Marker',marker,'LineWidth',1.5);
            end
            legend([{'GT'}, arrayfun(@(k)sprintf('C%d',k),1:c,'uni',false), {'Centers'}], ...
                   'Location','BestOutside');
            hold off;

            % build multiline title: heading + centers
            info = cell(c+1,1);
            info{1} = sprintf('%s (raw), Seed: %d', algMode, seed);
            for k=1:c
                cnt = sum(obj.labels==k);
                ctr = obj.centers(k,:);
                info{k+1} = sprintf('C%d: (%.2f, %.2f), n=%d', k, ctr(1), ctr(2), cnt);
            end
            title(info, 'Interpreter','none');

            %% Plot Pruned
            subplot(1,2,2);
            hold on;
            grid minor;

            %% Scatter GT
            scatter(gtMeans(:,1), gtMeans(:,2), 100, 'k', 's','LineWidth',1.5);

            for k=1:c
                scatter(X(obj.cleanLabels==k,1), X(obj.cleanLabels==k,2), 36, cols(k,:), 'filled');
                cnt = sum(obj.cleanLabels==k);
                marker = ternary(cnt>0, '*', 'x');
                scatter(obj.centers(k,1), obj.centers(k,2), 100, cols(k,:), ...
                        'Marker',marker,'LineWidth',1.5);
            end
            scatter(X(obj.outlierIdx,1), X(obj.outlierIdx,2), 36, [.5 .5 .5], 'x');

            legend([{'GT'}, arrayfun(@(k)sprintf('C%d',k),1:c,'uni',false), ...
                    {'Outliers','Centers'}], 'Location','BestOutside');
            hold off;

            % Cluster info: center and count per cluster
            info = cell(c+2,1);
            info{1} = sprintf('%s (pruned), Seed: %d', algMode, seed);
            for k=1:c
                cnt = sum(obj.cleanLabels==k);
                ctr = obj.centers(k,:);
                info{k+1} = sprintf('C%d: (%.2f, %.2f), n=%d', k, ctr(1), ctr(2), cnt);
            end
            title(info, 'Interpreter','none');

            plotbrowser('on');
        end
    end
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
%% ~~~~~~~~~~~~ Helper Functions ~~~~~~~~~~~~~~~ %%
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
function out = ternary(cond, a, b)
    if cond
        out = a; 
    else 
        out = b; 
    end
end
