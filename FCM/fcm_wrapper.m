function fcm_wrapper()
%DEMO_FCM_OUTLIERREMOVAL  Fuzzy C-Means clustering with per‐cluster outlier pruning
%
%   This demo:
%    1) Generates 3 Gaussian clusters in 2D
%    2) Runs vectorized FCM
%    3) Discards the top (1–x) fraction of points farthest from each center
%    4) Plots before and after pruning, with a legend that adapts
%       to the number of clusters.

    %% 1) Generate synthetic data
    rng(0);
    X1 = randn(100,2)*0.3 + [ 1, -1];
    X2 = randn(120,2)*0.4 + [-1,  2];
    X3 = randn( 80,2)*0.5 + [ 2,  1];
    X  = [X1; X2; X3];   % N×2

    %% 2) Run FCM
    c      = 3;                       % # clusters
    m      = 2.0;                     % fuzzifier
    metric = @(A,B) pdist2(A,B);      % Euclidean distance
    MaxIt  = 100;                     % max iterations
    tol    = 1e-5;                    % center‐shift tolerance
    [labels, centers] = fcm(c, X, m, metric, MaxIt, tol);

    %% 3) Discard outliers per cluster
    x = 0.90;  % keep the 90% closest points in each cluster
    [cleanLabels, isOutlier] = discardClusterOutliers(X, labels, centers, x);

    %% 4) Plot results
    colors = lines(c);
    figure('Position',[100 100 800 400]);

    %--- Before pruning ---
    subplot(1,2,1); hold on; axis equal; grid on;
    for k = 1:c
        idx = (labels==k);
        scatter(X(idx,1), X(idx,2), 36, colors(k,:), 'filled');
    end
    scatter(centers(:,1), centers(:,2), 100, 'k*','LineWidth',1.5);
    title('FCM Clustering (raw)');
    xlabel('X_1'); ylabel('X_2');
    legendEntries = [arrayfun(@(k)sprintf('Cluster %d',k), 1:c, 'Uni',false), {'Centers'}];
    legend(legendEntries,'Location','BestOutside');

    %--- After pruning ---
    subplot(1,2,2); hold on; axis equal; grid on;
    for k = 1:c
        idx = (cleanLabels==k);
        scatter(X(idx,1), X(idx,2), 36, colors(k,:), 'filled');
    end
    scatter(X(isOutlier,1), X(isOutlier,2), 36, [0.5 0.5 0.5], 'x');  % outliers
    scatter(centers(:,1), centers(:,2), 100, 'k*','LineWidth',1.5);
    title(sprintf('Pruned to %.0f%% per cluster, removed %d pts', x*100, sum(isOutlier)));
    xlabel('X_1'); ylabel('X_2');
    legendEntries = [arrayfun(@(k)sprintf('Cluster %d',k), 1:c, 'Uni',false), {'Outliers','Centers'}];
    legend(legendEntries,'Location','BestOutside');
end