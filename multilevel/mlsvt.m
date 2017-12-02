function [ Mapprox, UH, SH, VH, svp, logData ] = mlsvt( M, svdfn, svdParams,...
    threshold, iter, opts )
%% Multilevel Singular Value Thresholding
% M                         given matrix to be approximated
% svdfn                     function for computing svd
% svdParams                 cell array containing extra input parameters for @svd
% threshold                 function for soft/hard thresholding
% opts                      extra input parameters for restriction and
%                           prolongation
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
[ MCoarse, model, opts ] = restrictVarRPCA( M, varName, opts, mod(iter, 2)+1 );
%thresholdVal = thresholdVal / (model.operatorParams.A.prolongCoeff ^ (opts.depth - 1));
[ UH, SH, VH ] = svdfn( MCoarse, svdParams{:} );
coeff = (size(M, 2) / size(MCoarse, 2));
diagSH = diag(SH);

 %sqrt(max(opts.scale.(varName)) / 0.35^(opts.depth - 1));
svp = length(find(diagSH .* coeff > threshold));
diagSigmaH = diagSH(1:svp) - threshold ./ coeff;
diagSigmaH = diagSigmaH(diagSigmaH > 0);
svp = length(diagSigmaH);
SigmaH = diag(diagSigmaH);

varDims = [svp, size(M, 2)];
Vh = prolongateVar( VH(:, 1:svp)', model, opts, varName, varDims )';

Mapprox = UH(:, 1:svp) * SigmaH * Vh';

logData = [];


end

