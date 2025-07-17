classdef cHDBSCAN < cClusterer
    %cHDBSCAN  Hierarchical Density‐Based Spatial Clustering (pure MATLAB)
    %
    % Tunable Properties:
    %   minpts           – minimum # neighbors for core‐distance (default=5)
    %   minClusterSize   – drop any cluster smaller than this (default=5)
    %
    % After Fit(X):
    %   centers     – K×d medoids of each cluster
    %   labels      – N×1 cluster labels (0 = noise)
    %   cleanLabels – same as labels
    %   outlierIdx  – logical mask of labels==0

    properties
        minPts = 5;
        minClusterSize = 5;
    end

    methods
        %% ~~~~~~~~~~~~~~~ cHDBSCAN ~~~~~~~~~~~~~~~ %%
        function obj = cHDBSCAN(varargin)
            %NAME/VALUE constructor
            if mod(numel(varargin),2)~=0
                error('cHDBSCAN: need name/value pairs');
            end
            for k=1:2:numel(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k+1};
                else
                    error('cHDBSCAN: unknown property "%s"',varargin{k});
                end
            end
        end

        %% ~~~~~~~~~~~~~~~ Fit ~~~~~~~~~~~~~~~ %%
        function Fit(obj, X)

            %% hdbscan
            [labels, medoids] = obj.hdbscan(X);

            % drop too‐small clusters
            maxL  = max(labels);
            counts = histcounts(labels, 0.5:1:(maxL+0.5));
            small  = find(counts < obj.minClusterSize);
            for c = small
                labels(labels==c) = 0;
            end

            % store into base‐class properties
            obj.labels     = labels(:);
            obj.centers    = medoids;
            obj.cleanLabels = obj.labels;
            obj.outlierIdx  = (obj.labels==0);
        end

        %% ~~~~~~~~~~~~~~~ hdbscan ~~~~~~~~~~~~~~~ %%
        function [labels, medoids] = hdbscan(obj, X)
            %HDBSCAN  Hierarchical Density‐Based Spatial Clustering of Applications with Noise
            %
            %   [labels, medoids] = hdbscan(X, minPts)
            %
            %   Inputs:
            %     X       – N×d data matrix (N points in d dimensions)
            %     minPts  – minimum number of neighbors for core‐distance (and minimum cluster size)
            %
            %   Outputs:
            %     labels  – N×1 integer cluster labels (0 = noise)
            %     medoids – K×d coordinates of the medoid of each cluster (K clusters)
            %
            %   This is a pure‐MATLAB implementation:
            %     1) Compute core‐distances via the minPts‐th nearest neighbor.
            %     2) Build the mutual‐reachability graph.
            %     3) Construct the MST of that graph.
            %     4) Condense the hierarchy by removing edges from highest weight down.
            %     5) Compute cluster stability and extract the flat clustering.
            %     6) Compute medoids for each retained cluster.

            if nargin<2 || isempty(obj.minPts)
                error('hdbscan requires X (N×d) and minPts (scalar).');
            end
            [N, d] = size(X);
            if N < obj.minPts
                error('Number of points N must be >= minPts.');
            end

            %% 1) Pairwise distances and core distances
            % Requires Statistics Toolbox for pdist / squareform
            D = squareform( pdist(X) );            % N×N full distance matrix
            coreDist = zeros(N,1);
            for i = 1:N
                di = sort( D(i,:) );
                coreDist(i) = di( min(obj.minPts, numel(di)) );
            end

            %% 2) Mutual‐reachability matrix
            % mutual reachability: max{ coreDist(i), coreDist(j), D(i,j) }
            CR = max( coreDist * ones(1,N), ones(N,1) * coreDist' );
            mReach = max( CR, D );

            %% 3) Minimum‐spanning tree of mutual‐reachability graph
            G = graph( mReach );
            T = minspantree( G, 'Method','sparse' );
            edges = T.Edges;                       % table with EndNodes & Weight
            u = edges.EndNodes(:,1);
            v = edges.EndNodes(:,2);
            w = edges.Weight;

            % Sort edges by descending weight (highest distance first)
            [wSorted, idx] = sort( w, 'descend' );
            uS = u(idx);
            vS = v(idx);

            % Build symmetric adjacency for current cluster
            A = sparse( [u; v], [v; u], 1, N, N );

            %% 4) Hierarchical splitting
            % Each entry in `clusters` stores members + birth/death lambdas
            clusters = struct();
            clusters.members = 1:N;
            clusters.birth   = 0;
            clusters.death   = Inf;
            pointCluster = ones(N,1);

            for e = 1:numel(wSorted)
                i = uS(e); j = vS(e);
                cid = pointCluster(i);
                if pointCluster(j) ~= cid
                    continue;  % already split
                end
                lambda = 1 / wSorted(e);
                % record death of current cluster
                clusters(cid).death = lambda;
                % remove the edge
                A(i,j) = 0;  A(j,i) = 0;
                % find connected components of this cluster
                mem = clusters(cid).members;
                subA = A(mem, mem);
                comps = conncomp( graph(subA), 'Type','weak' );
                mask1 = (comps == comps(1));
                members1 = mem(mask1);
                members2 = mem(~mask1);
                % update parent cluster
                clusters(cid).members = members1;
                clusters(cid).birth   = lambda;
                % create new child cluster
                newID = numel(clusters) + 1;
                clusters(newID) = struct( ...
                    'members', members2, ...
                    'birth',   lambda, ...
                    'death',   Inf );
                % reassign points
                pointCluster(members1) = cid;
                pointCluster(members2) = newID;
            end

            %% 5) Compute stability and extract flat clustering
            K = numel(clusters);
            stability = zeros(K,1);
            for k = 1:K
                sz = numel(clusters(k).members);
                if sz >= obj.minPts
                    stability(k) = sz * (clusters(k).death - clusters(k).birth);
                end
            end
            [~, order] = sort(stability, 'descend');

            labels = zeros(N,1);
            nextLabel = 0;
            for id = order'
                mem = clusters(id).members;
                if numel(mem) < obj.minPts
                    continue;
                end
                unassigned = mem( labels(mem)==0 );
                if isempty(unassigned)
                    continue;
                end
                nextLabel = nextLabel + 1;
                labels(unassigned) = nextLabel;
            end

            %% 6) Compute medoids for each cluster
            uniqLabels = setdiff(unique(labels), 0)';
            Kf = numel(uniqLabels);
            medoids = zeros(Kf, d);
            for k = 1:Kf
                lbl = uniqLabels(k);
                pts = X(labels == lbl, :);
                Dsub = pdist2(pts, pts);
                [~, mi] = min( sum(Dsub,2) );
                medoids(k,:) = pts(mi,:);
            end
        end

    end
end
