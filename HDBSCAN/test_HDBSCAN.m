function test_HDBSCAN()
% TEST_HDBSCAN  Example script to exercise the HDBSCAN class
%
%   This script:
%     1) Generates 3 Gaussian clusters in 2D + some uniform noise
%     2) Instantiates HDBSCAN, sets parameters
%     3) Runs the full pipeline (fit_model → get_best_clusters → get_membership)
%     4) Plots condensed cluster tree and clustered data
%     5) Generates new test points and uses predict() to assign clusters
%     6) Visualizes predictions, highlighting outliers
    
    if false
        bla = HDBSCAN([randn(15,2); [5,5]+ randn(15,2)]);
        bla.run_hdbscan( 2,3,[])
        bla.plot_clusters();
    end

    rng(123);  % for reproducibility

    %% 1) Generate synthetic data
    N1 = 150; mu1 = [ 2,  2]; Sigma1 = 0.3*eye(2);
    N2 = 100; mu2 = [-2,  2]; Sigma2 = 0.5*eye(2);
    N3 = 120; mu3 = [ 0, -2]; Sigma3 = 0.4*eye(2);
    % Gaussian clusters
    X1 = mvnrnd(mu1, Sigma1, N1);
    X2 = mvnrnd(mu2, Sigma2, N2);
    X3 = mvnrnd(mu3, Sigma3, N3);
    % Uniform noise
    N0 = 50;
    X0 = (rand(N0,2)-0.5)*8;
    X  = [X1; X2; X3; X0];

    %% 2) Instantiate and configure HDBSCAN
    clusterer = HDBSCAN(X);
    % clusterer.minpts        = 10;   % core-dist neighbor count
    % clusterer.minclustsize  = 15;   % minimum size of a valid cluster
    % clusterer.outlierThresh = 0.8;  % loosen outlier assignment

    %% 3) Run full clustering pipeline
    % this calls fit_model, get_best_clusters, and get_membership internally
    % set the last argument to true to auto-plot results
    clusterer.run_hdbscan(clusterer.minpts, ...
                          clusterer.minclustsize, ...
                          1, ...          % min # clusters = 1
                          1, ...          % dEps = 1 (no skipping)
                          clusterer.outlierThresh, ...
                          true);          % plotResults = true

    %% 4) Overlay medoid / core points
    hold on;
    cores = clusterer.corePoints;
    % corePoints is a cell array: each cell contains indices of core
    % so pick representative point per cluster
    for k = 1:numel(cores)
        cp = cores{k};
        plot(X(cp,1), X(cp,2), 'kp', 'MarkerSize', 12, 'MarkerFaceColor','y');
    end
    hold off;

    %% 5) Predict on new data
    M = 200;
    % new points sampled in the same domain
    Y = (rand(M,2)-0.5)*8;
    [newLabels, newProb, outliers] = clusterer.predict(Y);

    %% 6) Plot predictions
    figure; hold on;
    % assign colors to predicted clusters
    cmap = lines(max(newLabels)); 
    for i = 1:M
        if newLabels(i)==0
            plot(Y(i,1), Y(i,2), 'x', 'Color',[0.5 0.5 0.5]);
        else
            plot(Y(i,1), Y(i,2), '.', 'Color', cmap(newLabels(i),:));
        end
    end
    title('Predicted cluster labels on new points');
    xlabel('Y(:,1)'); ylabel('Y(:,2)');
    legend({'Outliers','Cluster 1','Cluster 2','...'}, 'Location','BestOutside');
    hold off;

    %% 7) Simple assertions
    assert(~isempty(clusterer.labels), 'No labels were assigned!');
    nClusters = numel(unique(clusterer.labels(clusterer.labels>0)));
    fprintf('Detected %d clusters in training data.\n', nClusters);
    fprintf('Predicted %d/%d new points as outliers (threshold %.2f).\n', ...
            numel(outliers), M, clusterer.outlierThresh);

end
