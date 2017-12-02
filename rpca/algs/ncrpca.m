function [ L_t, S_t, iters, logData ] = ncrpca( M, true_r, opts )
% This matlab code implements Non-convex Robust PCA (NcRPCA)
% Input:
% M = given low rank+sparse matrix to be decomposed
% true_r = maximum rank of the low rank rank component
% EPS (optional) = convergence threshold for ||M-(L_t+S_t)||_F; default is 1e-3
% MAX_ITER (optional) = maximum iterations for NcRPCA; default is 51
% EPS_S (optional) = threshold for removing small entries in the sparse component; default is 1e-3
% incoh (optional) = incoherence of the low rank component; default is 1
% TOL (optional) = tolerance for relative error in ||M-(L_t+S_t)||_F in consecutive iterations; default is 1e-1
% Output:
% M_t = thresholded M at each iteration
% L_t = rank-k approximation of M_t
% S_t = sparse component, computed as M-M_t
% iters = number of iteration of NcRPCA
% frob_err = ||M-(L_t+S_t)||_F at each iteration

startTime = cputime;

[~, n] = size(M);
t = 1;
idx = [];
thresh_const = opts.incoh; % threshold constant: can be tuned depending on incoherence
thresh_red = 0.9; % parameter to reduce the threshold constant: can be tuned
r_hat = 1; % initial rank for stagewise algorithm
L_t = zeros(size(M));
Sig_t = lansvd(M,1,'L');
D_t = M-L_t;
thresh = thresh_const*Sig_t/sqrt(n);
idx = unique([find(abs(D_t) > thresh); idx]);
S_t = zeros(size(M));
S_t(idx) = D_t(idx); % initial thresholding
if max(idx(:))==0
    idx = [];
end
normM = norm(M, 'fro');

logData.frob_err = zeros(opts.maxIter, 1);
logData.frob_err(1) = inf;
logData.errL = zeros(opts.maxIter, 1);
logData.errS = zeros(opts.maxIter, 1);
logData.timeSteps = zeros(opts.maxIter, 1);
logData.timeSteps(1) = startTime;
timeStep = 1;

while logData.frob_err(t) / normM >= opts.tolerance && t<opts.maxIter % convergence check
    t = t+1;
    [U_t, Sig_t, V_t] = lansvd(M-S_t, r_hat+1, 'L');
    %[U_t, Sig_t, V_t] = svds(M-S_t, r_hat+1); % use this if not using propack
    L_t= U_t(:,1:r_hat)*Sig_t(1:r_hat,1:r_hat)*V_t(:,1:r_hat)';
    D_t = M-L_t;
    thresh = (thresh_const/sqrt(n))*Sig_t(r_hat+1, r_hat+1); % use n instead of sqrt(n) for thresholding less aggresively
    idx = unique([find(abs(D_t) > thresh); idx]);
    S_t(idx) = D_t(idx);
    logData.frob_err(t) = norm(M-(L_t+S_t), 'fro');
    if ~mod(t, opts.logInterval)
        disp(['Iter no. ' num2str(t) '; err = ' num2str(logData.frob_err(t) / normM)]);
    end
    if ((logData.frob_err(t-1)-logData.frob_err(t))/logData.frob_err(t-1) <= opts.TOL) && r_hat<true_r
%         r_hat = r_hat+1; % use this for incrementally updating rank by 1
        sig_t = lansvd(M-S_t, true_r, 'L'); % svd function from propack
        ratio_sig = sig_t(r_hat+1:end)./[sig_t(r_hat+2:end); sig_t(end)];
        [~, mx_idx] = max(ratio_sig);
        r_hat = r_hat+mx_idx; % update rank for the next stage
    elseif ((logData.frob_err(t-1)-logData.frob_err(t))/logData.frob_err(t-1) <= opts.TOL) && r_hat==true_r
        thresh_const = thresh_const*thresh_red; % tune threshold
    end
    
    currTime = cputime;
    if currTime - logData.timeSteps(timeStep) >= opts.timeStepLen
        logData.timeSteps(timeStep) = currTime - startTime;
        if isfield(opts, 'L0') && ~isempty(opts.L0)
            logData.errL(timeStep) = norm(L_t - opts.L0, 'fro') / norm(opts.L0, 'fro');
        end
        if isfield(opts, 'S0') && ~isempty(opts.S0)
            logData.errS(timeStep) = norm(S_t - opts.S0, 'fro') / norm(opts.S0, 'fro');
        end
        logData.timeSteps(timeStep + 1) = logData.timeSteps(timeStep);
        timeStep = timeStep + 1;
    end
end
S_t(abs(S_t)<opts.EPS_S) = 0; % threshold to remove small errors and obtain sparse component
iters = t; % no. of iters. of ncrpca = length(frob_err)-1

disp([' == Finishing with: Iter no. ' num2str(t) '; err = ' num2str(logData.frob_err(t) / normM)]);

logData.frob_err = logData.frob_err(1 : t);
logData.errL = logData.errL(1 : timeStep - 1);
logData.errS = logData.errS(1 : timeStep - 1);
logData.timeSteps = logData.timeSteps(1 : timeStep - 1);

end
