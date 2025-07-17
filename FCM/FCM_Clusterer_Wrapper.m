function FCM_Clusterer_Wrapper()
% Generate demo data
algMode = "DBSCAN";
seedVec = [1:10];

gtMeans = [...
    1,  -1;
    -1, 2;
    2,  1;
    3,  4;
    9,  3;
    ];

for seed = seedVec
    rng(seed);
    X1 = randn(5,2)*0.3 + gtMeans(1,:);
    X2 = randn(4,2)*0.4 + gtMeans(2,:);
    X3 = randn(3,2)*0.5 + gtMeans(3,:);
    X4 = randn(2,2)*0.5 + gtMeans(4,:);
    X5 = randn(1,2)*0.5 + gtMeans(5,:);
    X  = [X1; X2; X3; X4; X5];

    switch algMode
        case "DBSCAN"
            oDBSCAN = cDBSCAN( ...
                'minClusterSize',   1, ...
                'maxDist',          1);

            %% Fit
            oDBSCAN.Fit(X);

            %% Plot
            oDBSCAN.Plot_results(X, gtMeans, seed);
        case "FCM"

            %% Create clusterer with automatic c‐selection
            oFCM = cFCM( ...
                'autoSelect',           true, ...
                'clusterRange',         1:8, ...
                'minClusterSize',       3, ...
                'validityIndex',        'XB', ...
                'fuzzifier',            5, ...
                'shapePenaltyWeight',   0.2, ...
                'outlierMethod',        'iqr', ...
                'outlierIQRFactor',     1.2);

            % profile on

            %% Fit
            oFCM.Fit(X);

            %% Plot
            oFCM.Plot_results(X, gtMeans, seed);
    end
end

% profile viewer

end
