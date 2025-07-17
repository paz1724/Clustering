classdef cDBSCAN < cClusterer
    
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~ Properties ~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    properties
        minClusterSize     = 3;
        maxDist            = 2;
        numClusters        = []

    end

    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~ Methods ~~~~~~~~~~~~~~~~~~~ %%
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
    methods
        %% ~~~~~~~~~~~~~~~ FCMClusterer ~~~~~~~~~~~~~~~ %%
        function obj = cDBSCAN(varargin)
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

            obj.labels  = dbscan(X,obj.maxDist,obj.minClusterSize);
            classes     = unique(obj.labels(obj.labels>0));
            N           = numel(classes);
            obj.centers = zeros(N,2);
            for i = 1:numel(classes)
                mask             = obj.labels==i;
                obj.centers(i,:) = mean(X(mask,:),1);
            end
            obj.cleanLabels = obj.labels;
            obj.outlierIdx  = find(obj.labels==-1);
            obj.numClusters = numel(classes);

        end

    end
end
