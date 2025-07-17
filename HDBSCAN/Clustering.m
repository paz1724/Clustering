% Example data
rng(0);
X = [randn(100,2)*0.5 + [2,2];
     randn( 80,2)*0.3 + [-2,2];
     randn( 60,2)*0.4 + [0,-2]];

minPts = 10;
[labels, medoids] = hdbscan(X, minPts);

% Plot
figure; hold on;
gscatter(X(:,1), X(:,2), labels);
plot(medoids(:,1), medoids(:,2), 'kp', 'MarkerSize',15,'MarkerFaceColor','y');
title('HDBSCAN (no toolboxes)');
