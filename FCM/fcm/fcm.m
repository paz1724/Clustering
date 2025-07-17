function [prediction, v] = fcm(c, X, m, metric, MaxIt, tol)
%FCM  Vectorized Fuzzy C-Means clustering
%   [prediction, v] = fcm(c, X, m, metric, MaxIt, tol)

    [N, D] = size(X);
    p      = 2/(m-1);

    % initialize membership U (c×N)
    U = rand(c,N);
    U = bsxfun(@rdivide, U, sum(U,1));

    % initial centers
    Um   = U.^m;                       % c×N
    v    = (Um * X) ./ (sum(Um,2)*ones(1,D));  % c×D
    vOld = v;
    delta = inf;
    it    = 0;

    while it<MaxIt && delta>tol
        % distance matrix (N×c)
        Dmat = metric(X, v);

        % update U
        DP    = Dmat .^ (-p);          % N×c
        denom = sum(DP,2);             % N×1
        U     = bsxfun(@rdivide, DP', denom');  % c×N

        % update centers
        Um   = U.^m;                   % c×N
        v    = (Um * X) ./ (sum(Um,2)*ones(1,D));  % c×D

        % check convergence
        delta = max(abs(v(:)-vOld(:)));
        vOld  = v;
        it    = it + 1;
    end

    % hard labels
    [~, prediction] = max(U, [], 1);
    prediction = prediction(:);
end