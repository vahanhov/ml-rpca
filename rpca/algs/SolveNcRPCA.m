function [ xStar, FStar, nIter, logData, model, opts ] = SolveNcRPCA(...
    model, params, opts )
%% Wrapper function for ncrpca
% 
%   Author: Vahan Hovhannisyan, 2017.


opts = setDefaultOption(opts, 'EPS_S', 1e-3); 
opts = setDefaultOption(opts, 'incoh', 1);
opts = setDefaultOption(opts, 'TOL', 1e-1);

[ xStar.L, xStar.S, nIter, logData ] = ncrpca( params.D, params.rank, opts ); %, EPS_S, incoh, TOL);

FStar = nan;


end

