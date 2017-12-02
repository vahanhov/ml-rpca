function [ Mapprox, UH, SH, VH ] = mlRankApprox( M, svdfn, svdParams, opts )
%% Calculates a low rank approximation for the given matrix A
% A                 given matrix to be approximated
% svd               function for computing svd
% svdParams         cell array containing extra input parameters for @svd
% opts              cell array containing extra input parameters for
% restriction and prolongation
%
%   Author: Vahan Hovhannisyan, 2017.

if isnumeric(svdParams{1})
    oldDepth = opts.depth;
    opts.depth = min(opts.depth, floor(log2(size(M, 2) / svdParams{1}) + 1));
    if oldDepth ~= opts.depth
        disp([' ==== Depth changed to ' num2str(opts.depth) ' from ' num2str(oldDepth) ' ==== ']);
    end
end

varName = 'V';
[ MCoarse, model, opts ] = restrictVarRPCA( M, varName, opts );
%thresholdVal = thresholdVal / (model.operatorParams.A.prolongCoeff ^ (opts.depth - 1));
[ UH, SH, VH ] = svdfn( MCoarse, svdParams{:} );

threshold = svdParams{1} - 1;
varDims = [threshold, size(M, 2)];
Vh = prolongateVar( VH(:, 1:threshold)', model, opts, varName, varDims )';

Mapprox = UH(:, 1:threshold) * SH(1:threshold, 1:threshold) * Vh(:, 1:threshold)';


end

