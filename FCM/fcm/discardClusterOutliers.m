function [cleanLabels, outlierIdx] = discardClusterOutliers(X, labels, centers, x)
    labels      = labels(:);
    N           = size(X,1);
    K           = size(centers,1);
    cleanLabels = labels;
    outlierIdx  = false(N,1);

    for k = 1:K
        idx = find(labels==k);
        if isempty(idx), continue; end

        % distances to center k
        diffs = X(idx,:) - centers(k,:);
        dists = sqrt(sum(diffs.^2,2));

        % threshold at x-fraction
        thr = prctile(dists, x*100);

        % mark outliers
        bad = dists > thr;
        outlierIdx(idx(bad))  = true;
        cleanLabels(idx(bad)) = 0;
    end
end
