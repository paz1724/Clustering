classdef cHDBSCAN < cClusterer
    %cDBSCAN  Hierarchical Density‐Based Spatial Clustering (HDBSCAN) in pure MATLAB
    %
    % Tunable Properties:
    %   minPts           – minimum # neighbors for core‐distance and cluster size
    %   minClusterSize   – drop any cluster smaller than this (default = minPts)
    %
    % After Fit(X):
    %   centers    – K×d medoids of each surviving cluster
    %   labels     – N×1 integer cluster labels (0 = noise)
    %   cleanLabels– same as labels
    %   outlierIdx – logical mask (labels==0)
    %   numClusters– number of non‐noise clusters

    properties
        minPts           = 5
        minClusterSize   = []   % if empty, will default to minPts
    end

    methods
        %% Constructor
        function obj = cHDBSCAN(varargin)
            % Allow name/value pairs
            if mod(numel(varargin),2)~=0
                error('cDBSCAN: need name/value pairs');
            end
            for k=1:2:numel(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k+1};
                else
                    error('cDBSCAN: unknown property "%s"',varargin{k});
                end
            end
            if isempty(obj.minClusterSize)
                obj.minClusterSize = obj.minPts;
            end
        end

        %% Fit
        function Fit(obj, X)
            %FIT  Run HDBSCAN on X (N×d)
            [N, d] = size(X);
            mpts   = obj.minPts;

            % 1) Distance matrix & core distances (exclude self)
            D = squareform( pdist(X) );
            coreDist = zeros(N,1);
            for i=1:N
                others = [1:i-1, i+1:N];
                di     = sort( D(i,others) );
                coreDist(i) = di( min(mpts, numel(di)) );
            end

            % 2) Mutual‐reachability
            CR    = max( coreDist*ones(1,N), ones(N,1)*coreDist' );
            mReach= max( CR, D );

            % 3) MST
            G  = graph(mReach);
            T  = minspantree(G,'Method','sparse');
            E  = T.Edges;
            [wSorted, idx] = sort(E.Weight,'descend');
            uS = E.EndNodes(idx,1);
            vS = E.EndNodes(idx,2);

            % adjacency for splits
            A = sparse([uS;vS],[vS;uS],1,N,N);

            % 4) Initialize hierarchy
            clusters      = struct('members',{{1:N}},'birth',0,'death',Inf);
            pointCluster  = ones(N,1);

            % 5) Split by removing edges in descending weight
            for e=1:numel(wSorted)
                i = uS(e); j = vS(e);
                cid = pointCluster(i);
                if pointCluster(j)~=cid, continue; end
                lambda = 1/wSorted(e);
                clusters(cid).death = lambda;

                % remove edge
                A(i,j)=0; A(j,i)=0;

                % connected comps within this cluster
                mem  = clusters(cid).members;
                subA = A(mem,mem);
                comp = conncomp(graph(subA),'Type','weak');
                mask1= (comp==comp(1));
                m1   = mem(mask1);
                m2   = mem(~mask1);

                % update parent
                clusters(cid).members = m1;
                clusters(cid).birth   = lambda;
                % new child
                newID = numel(clusters)+1;
                clusters(newID) = struct(...
                    'members',{m2},...
                    'birth',   lambda,...
                    'death',   Inf);

                % reassign points
                pointCluster(m1) = cid;
                pointCluster(m2) = newID;
            end

            % 6) Compute stability & extract flat clustering
            K = numel(clusters);
            stab = zeros(K,1);
            for k=1:K
                sz = numel(clusters(k).members);
                if sz >= mpts
                    stab(k) = sz * (clusters(k).death - clusters(k).birth);
                end
            end
            [~, order] = sort(stab,'descend');

            labels = zeros(N,1);
            nextL  = 0;
            for id = order'
                mem = clusters(id).members;
                if numel(mem) < mpts, continue; end
                unassigned = mem(labels(mem)==0);
                if isempty(unassigned), continue; end
                nextL = nextL + 1;
                labels(unassigned) = nextL;
            end

            % 7) Drop small clusters
            good = false(max(labels),1);
            for k=1:max(labels)
                if sum(labels==k) >= obj.minClusterSize
                    good(k) = true;
                end
            end
            for k=1:max(labels)
                if ~good(k)
                    labels(labels==k) = 0;
                end
            end

            % 8) Compute medoids of surviving clusters
            uniqL = setdiff(unique(labels),0);
            Kf    = numel(uniqL);
            medoids = zeros(Kf,d);
            for ii=1:Kf
                pts = X(labels==uniqL(ii),:);
                Dsub= pdist2(pts,pts);
                [~,mi] = min(sum(Dsub,2));
                medoids(ii,:) = pts(mi,:);
            end

            % Store results
            obj.labels      = labels;
            obj.cleanLabels = labels;
            obj.outlierIdx  = labels==0;
            obj.centers     = medoids;
            obj.numClusters = Kf;
        end
    end
end
