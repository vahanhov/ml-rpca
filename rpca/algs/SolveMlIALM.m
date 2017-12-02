function [ xStar, FStar, nIter, logData, model, opts ] = SolveMlIALM(...
    model, params, opts )
%% Wrapper function for inexact_alm_rpca_ml
% 
%   Author: Vahan Hovhannisyan, 2017.

disp(['-- using up to ' num2str(opts.depth) ' levels for multi-level SVDs --']);

[ xStar.L, xStar.S, nIter, logData ] = inexact_alm_rpca_ml( params.D, params.nu, opts ); %, opts.tolerance, opts.maxIter );

FStar = nan;


end
