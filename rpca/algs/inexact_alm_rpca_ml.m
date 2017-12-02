function [ A_hat, E_hat, iter, logData ] = inexact_alm_rpca_ml( D, lambda0, opts )
%% Based on the inexact_alm_rpca function, with the singular value thresholding
%   replaced by a multilevel approximation
%
%   Author: Vahan Hovhannisyan, 2017.

startTime = cputime;

[m, n] = size(D);

% initialize
lambda = lambda0;% / opts.lambdaCoeff;
if isfield(opts, 'lambdaMlCoeff')
    lambda = lambda * opts.lambdaMlCoeff;
end

Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);

if isfield(opts, 'mu') && ~isempty(opts.mu)
    mu = opts.mu;
else
    %[ YCoarse, ~, ~ ] = restrictVarRPCA( Y, opts.coarseVarNames{1}, opts );
    mu = 5 / norm_two; %lansvd(YCoarse, 1, 'L');%_coarse; %norm_two;%_coarse; % this one can be tuned
end
mu_bar = 1e7 * mu;
if isfield(opts, 'rho') && ~isempty(opts.rho)
    rho = opts.rho;
else
    rho = 1.5;          % this one can be tuned
end
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
sv = 10;

opts.depth = min(opts.depth, floor(log(n / sv) + 1));
    
logData.rank_iters = -ones(opts.maxIter, 1);
logData.frob_err = -ones(opts.maxIter, 1);
logData.S1_iters = -ones(opts.maxIter, 1);
logData.errL = -ones(opts.maxIter, 1);
logData.errS = -ones(opts.maxIter, 1);
logData.timeSteps = -ones(opts.maxIter, 1);
logData.timeSteps(1) = startTime;
mlsvt_max = 5;
logData.mlsvt = -ones(opts.maxIter, mlsvt_max);
timeStep = 1;

while ~converged
    iter = iter + 1;
    
%     temp_T = D - A_hat + (1/mu)*Y;
%     E_hat = max(temp_T - lambda/mu, 0);
%     E_hat = E_hat + min(temp_T + lambda/mu, 0);
    
    %% Update E   
%     if doCoarse
%         sparseThreshold = lambda/mu;
%     else
%         sparseThreshold = lambda/mu;
%     end
    E_hat = soft(D - A_hat + (1/mu)*Y, lambda/mu); % * S1); %sparseThreshold);

    %% Singulat value thresholding
    logData.rank_iters(iter) = sv;
    
    if choosvd(n, sv) == 1
        svdfn = @lansvd;
        svdParams = {sv, 'L'};
    else
        svdfn = @svd;
        svdParams = {'econ'};
    end
    
    [ A_hat, U, S, V, svp, logData_mlsvt ] = mlsvt( D - E_hat + (1/mu)*Y,...
        svdfn, svdParams, 1/mu, iter, opts );
    
    if ~isempty(logData_mlsvt)
        mlsvt_len = min(length(logData_mlsvt), mlsvt_max);
        logData.mlsvt(iter, 1:mlsvt_len) = logData_mlsvt(:, 1:mlsvt_len);
    end
    coarseInfo = '';
    logData.svp(iter) = svp;
    logData.sv(iter) = sv;
    % svp is the number of non-zero singular values after thresholding
    % sv is the number of singular values to compute
    
    %% If the coarse level was not useful, do a fine iteration
    %   This results in starting with fine iterations and finishing with coarse
    if false% iter > 10 && all(logData.svp(end-10 : end) == svp) %false%svp < 1 %sv / 1.5; %mu <= sqrt(n / size(Ycoarse, 2));
        if choosvd(n, sv) == 1
            [U, S, V] = lansvd(D - E_hat + (1/mu)*Y, sv, 'L');
        else
            [U, S, V] = svd(D - E_hat + (1/mu)*Y, 'econ');
        end
        diagS = diag(S);
        svp = length(find(diagS > 1/mu));

        A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';
        coarseInfo = '; fine iteration';
    end
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    logData.S1_iters(iter) = S(1, 1);
    
    total_svd = total_svd + 1;
    %% Update Y and mu
    Z = D - A_hat - E_hat;
    Y = Y + mu*Z;
    mu = min(rho * mu, mu_bar);
    %doCoarse = mu ~= mu_bar;
    
    %% stop Criterion
    logData.frob_err(iter) = norm(Z, 'fro') / d_norm;
    stopCriterion = logData.frob_err(iter);
    if opts.maxTime <= 0
        if stopCriterion < opts.tolerance
            converged = true;
        end
    end
    
    currTime = cputime;
    if currTime - logData.timeSteps(timeStep) >= opts.timeStepLen
        logData.timeSteps(timeStep) = currTime - startTime;
        if isfield(opts, 'L0') && ~isempty(opts.L0)
            logData.errL(timeStep) = norm(A_hat(:, opts.GTidx) - opts.L0, 'fro') / norm(opts.L0, 'fro');
        end
        if isfield(opts, 'S0') && ~isempty(opts.L0)
            logData.errS(timeStep) = norm(E_hat(:, opts.GTidx) - opts.S0, 'fro') / norm(opts.S0, 'fro');
        end
        logData.timeSteps(timeStep + 1) = logData.timeSteps(timeStep);
        timeStep = timeStep + 1;
    end
    
    if mod( total_svd, 10) == 0 %|| ~isempty(coarseInfo)
        disp(['#svd ' num2str(total_svd) ...%' r(A) ' num2str(rank(A_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)) / numel(E_hat))...
            ' stopCriterion ' num2str(stopCriterion)...
            ' mu = ' num2str(mu)...
            coarseInfo]);
    end    
    
    if ~converged && iter >= opts.maxIter
        disp('Maximum iterations reached');
        converged = 1;
    end
    
    if opts.maxTime > 0 && cputime - startTime >= opts.maxTime
        disp('Maximum time reached');
        converged = 1;
    end
end

% E_hat = abs(E_hat);
E_hat(abs(E_hat) < opts.epsRound) = 0; % threshold to remove small errors and obtain sparse component

logData.frob_err = logData.frob_err(1 : iter);
logData.rank_iters = logData.rank_iters(1 : iter);
logData.S1_iters = logData.S1_iters(1 : iter);
logData.errL = logData.errL(1 : timeStep - 1);
logData.errS = logData.errS(1 : timeStep - 1);
logData.timeSteps = logData.timeSteps(1 : timeStep - 1);

logData.mlsvt = logData.mlsvt(1 : iter, :);

end