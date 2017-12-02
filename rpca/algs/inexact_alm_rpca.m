function [ A_hat, E_hat, iter, logData ] = inexact_alm_rpca( D, lambda, opts )
% Oct 2009
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust PCA.
%
% D - m x n matrix of observations/data (required input)
%
% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
% Initialize A,E,Y,u
% while ~converged 
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_* + lambda * |E|_1 + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end
%
% Minming Chen, October 2009. Questions? v-minmch@microsoft.com ; 
% Arvind Ganesh (abalasu2@illinois.edu)
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

%addpath PROPACK;

startTime = cputime;

[m, n] = size(D);

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
if isfield(opts, 'mu') && ~isempty(opts.mu)
    mu = opts.mu;
elseif isfield(opts, 'muCoeff') && ~isempty(opts.muCoeff)
    mu = opts.muCoeff / norm_two;%_coarse; %norm_two;%_coarse; % this one can be tuned
else
    mu = 1 / norm_two;
end
mu_bar = mu * 1e7;
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

logData.rank_iters = -ones(opts.maxIter, 1);
logData.frob_err = -ones(opts.maxIter, 1);
logData.S1_iters = -ones(opts.maxIter, 1);
logData.errL = -ones(opts.maxIter, 1);
logData.errS = -ones(opts.maxIter, 1);
logData.timeSteps = zeros(opts.maxIter, 1);
logData.timeSteps(1) = startTime;
timeStep = 1;

while ~converged
    iter = iter + 1;
    
    temp_T = D - A_hat + (1/mu)*Y;
    E_hat = max(temp_T - lambda/mu, 0);
    E_hat = E_hat+min(temp_T + lambda/mu, 0);
    %sv
    logData.rank_iters(iter) = sv;
    if choosvd(n, sv) == 1
        [U, S, V] = lansvd(D - E_hat + (1/mu)*Y, sv, 'L');
    else
        [U, S, V] = svd(D - E_hat + (1/mu)*Y, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    
    logData.svp(iter) = svp;
    logData.sv(iter) = sv;
    
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    logData.S1_iters(iter) = diagS(1);
    
    A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';    
    total_svd = total_svd + 1;
    
    Z = D - A_hat - E_hat;
    
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
        
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
        if isfield(opts, 'S0') && ~isempty(opts.S0)
            logData.errS(timeStep) = norm(E_hat(:, opts.GTidx) - opts.S0, 'fro') / norm(opts.S0, 'fro');
        end
        logData.timeSteps(timeStep + 1) = logData.timeSteps(timeStep);
        timeStep = timeStep + 1;
    end
    
    if mod( total_svd, 10 ) == 0
        disp(['#svd ' num2str(total_svd) ...%' r(A) ' num2str(rank(A_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)) / numel(E_hat))...
            ' stopCriterion ' num2str(stopCriterion)...
            ' mu = ' num2str(mu)]);
    end    
    
    if ~converged && iter >= opts.maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
    
    if opts.maxTime > 0 && cputime - startTime >= opts.maxTime
        disp('Maximum time reached');
        converged = 1;
    end
end

% E_hat = abs(E_hat);
% E_hat(abs(E_hat) < opts.epsRound) = 0; % threshold to remove small errors and obtain sparse component

logData.frob_err = logData.frob_err(1 : iter);
logData.rank_iters = logData.rank_iters(1 : iter);
logData.S1_iters = logData.errS(1 : iter);
logData.errL = logData.errL(1 : timeStep - 1);
logData.errS = logData.errS(1 : timeStep - 1);
logData.timeSteps = logData.timeSteps(1 : timeStep - 1);

end
