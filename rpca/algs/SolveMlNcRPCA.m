function [ xStar, FStar, nIter, logData, model, opts ] = SolveMlNcRPCA(...
    model, params, opts )
%% Wrapper function for ncrpca_ml
% 
%   Author: Vahan Hovhannisyan, 2017.

disp(['-- using up to ' num2str(opts.depth) ' levels for multi-level SVDs --']);

opts = setDefaultOption(opts, 'EPS_S', 1e-3); 
opts = setDefaultOption(opts, 'incoh', 1);
opts = setDefaultOption(opts, 'TOL', 1e-1);
if isfield(opts, 'rank')
    params.rank = opts.rank;
    disp([' == Changing the rank to ' num2str(opts.rank)]);
end

[ xStar.L, xStar.S, nIter, logData ] = ncrpca_ml( params.D, params.rank, opts ); %, EPS_S, incoh, TOL);

FStar = nan;


end

