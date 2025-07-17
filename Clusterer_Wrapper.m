function Clusterer_Wrapper()

%% Generate demo data
algModeList = ["MeanShift", "DBSCAN", "FCM"];    % ["MeanShift", "HDBSCAN", "DBSCAN", "FCM"]
seedVec = 6;
gtMeans = [
    1,  -1;
    -1,  2;
    2,   1;
    3,   4;
    9,   3;
    ];
gtNum = [...
    5; ...
    4; ...
    3; ...
    3; ...
    1; ...
    ];

gtStd = [...
    0.3; ...
    0.4; ...
    0.5; ...
    0.5; ...
    0.5; ...
    ];

for algMode = algModeList
    %% Create Clustering Object
    oCluster = Create_clustering_object(algMode);

    %% sPlot
    sPlot.gtMeans   = gtMeans;
    sPlot.algMode   = algMode;

    %% seed Loop
    for seed = seedVec
        rng(seed);

        %% Generate Input Data
        X = [];
        for i = 1:numel(gtNum)
            Xi = randn(gtNum(i),2)*gtStd(i) + gtMeans(i,:);
            X = vertcat(X, Xi);
        end

        %% Fit
        oCluster.Fit(X);

        %% Plot_results
        sPlot.seed      = seed;
        oCluster.Plot_results(X, sPlot);

    end

end

% profile viewer

end

%% ~~~~~~~~~~~ Create_clustering_object ~~~~~~~~~~~ %%
function [oCluster] = Create_clustering_object(algMode)
switch algMode
    case "DBSCAN"
        oCluster = cDBSCAN( ...
            'minClusterSize', 1, ...
            'maxDist',        1 );

    case "HDBSCAN"
        oCluster = cHDBSCAN( ...
            'minPts',          3, ...
            'minClusterSize',  3 );

    case "FCM"
        oCluster = cFCM( ...
            'autoSelect',        true, ...
            'clusterRange',      1:8, ...
            'minClusterSize',    3, ...
            'validityIndex',     'XB', ...
            'fuzzifier',         2, ...
            'shapePenaltyWeight',0, ...
            'outlierMethod',     'iqr', ...
            'outlierIQRFactor',  1.2 );

    case "MeanShift"
        oCluster = cMeanShift( ...
            'bandWidth',      1.0, ...
            'plotFlag',       false, ...
            'minClusterSize', 2 );
end

end
