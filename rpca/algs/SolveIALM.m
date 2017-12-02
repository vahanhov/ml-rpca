function [ xStar, FStar, nIter, logData, model, opts ] = SolveIALM(...
    model, params, opts )
%% Wrapper function for inexact_alm_rpca
% 
%   Author: Vahan Hovhannisyan, 2017.

[ xStar.L, xStar.S, nIter, logData ] = inexact_alm_rpca( params.D, params.nu, opts );

FStar = nan;


end
