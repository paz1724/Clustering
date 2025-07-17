classdef cMeanShift < cClusterer
    %cMEANSHIFT  Mean‐Shift clustering with flat kernel
    %
    % Tunable Properties:
    %   bandWidth        - kernel radius
    %   plotFlag         - true to show the per-iteration plot
    %   minClusterSize   - drop any cluster smaller than this
    %
    % After Fit(X):
    %   centers          - (c×d) cluster centers
    %   labels           - (N×1) hard labels (0=noise)
    %   cluster2dataCell - (c×1) cells of point‐indices

    properties
        bandWidth      = 1
        plotFlag       = false
        minClusterSize = 1
    end

    properties (SetAccess=protected)
        cluster2dataCell
    end

    methods
        %% ~~~~~~~~~~~~~~~ cMeanShift ~~~~~~~~~~~~~~~ %%
        function obj = cMeanShift(varargin)
            %NAME/VAL constructor
            if mod(numel(varargin),2)~=0
                error('cMeanShift: need name/value pairs');
            end
            for k=1:2:numel(varargin)
                if isprop(obj,varargin{k})
                    obj.(varargin{k}) = varargin{k+1};
                else
                    error('cMeanShift: unknown property "%s"',varargin{k});
                end
            end
        end

        %% ~~~~~~~~~~~~~~~ Fit ~~~~~~~~~~~~~~~ %%
        function Fit(obj, X)
            %FIT  Run MeanShiftCluster (expects d×N)
            [cc, lbl, c2d] = MeanShiftCluster( X', obj.bandWidth, obj.plotFlag );
            obj.centers = cc';              % c×d
            obj.labels  = lbl(:);           % N×1
            obj.cluster2dataCell = c2d;     % c×1
            obj.cleanLabels = obj.labels;

            obj.outlierIdx  = find(obj.labels==-1);

            %% drop too‐small clusters
            counts = cellfun(@numel, c2d);
            small  = find(counts < obj.minClusterSize);
            if ~isempty(small)
                mask = ismember(obj.labels, small);
                obj.cleanLabels(mask) = 0;
            end

            
        end

    end
end
