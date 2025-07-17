function hdbscan_wrapper()
%TEST_HDBSCAN  Test script for the hdbscan function

    %% 1) Generate synthetic data
    rng(42);
    N1 = 100; mu1=[2,2];  S1=0.3*eye(2);
    N2 =  80; mu2=[-2,2]; S2=0.4*eye(2);
    N3 =  60; mu3=[0,-2]; S3=0.5*eye(2);
    X1 = mvnrnd(mu1,S1,N1);
    X2 = mvnrnd(mu2,S2,N2);
    X3 = mvnrnd(mu3,S3,N3);
    N0 = 50; X0 = (rand(N0,2)-0.5)*8;
    X  = [X1; X2; X3; X0];

    %% 2) Run HDBSCAN
    minPts = 10;
    [labels, medoids] = hdbscan(X, minPts);

    %% 3) Plot results
    figure; hold on;
    maxL = max(labels);
    cmap = lines(maxL);
    % plot noise
    scatter(X(labels==0,1), X(labels==0,2), 36, 'k', 'x');
    % plot clusters
    for c = 1:maxL
        scatter(X(labels==c,1), X(labels==c,2), 36, cmap(c,:), 'filled');
    end
    % overlay medoids
    scatter(medoids(:,1), medoids(:,2), 100, 'kp', 'filled','MarkerEdgeColor','w');
    title('HDBSCAN Clustering Test');
    legend(['Noise', arrayfun(@(c) sprintf('Cluster %d',c),1:maxL,'uni',false), 'Medoids'], ...
           'Location','BestOutside');
    hold off;

    %% 4) Assertions
    assert(numel(labels)==size(X,1), 'Labels size mismatch.');
    assert(size(medoids,1)==maxL,   'Medoid count mismatch.');
    assert(size(medoids,2)==size(X,2),'Medoid dimension mismatch.');
    assert(maxL>=1,                  'No clusters detected.');
    fprintf('test_hdbscan: PASSED — %d clusters found (excluding noise).\n', maxL);
end