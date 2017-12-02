function [ Acoarse, model, opts ] = restrictVarRPCA_left( A, varName, varargin )
%% Creates a coarse matrix for the RPCA problems 
%   by applying the restriction operator from left
%
%   also returns model and opts  structs that should be passed to the
%   prolongation function later
%
%   Author: Vahan Hovhannisyan, 2017.


if ~isempty(varargin)
    opts = varargin{1};
else
    opts = [];
end

% if ~isfield(opts, 'depth')
%     opts.depth = floor(log2(size(A, 2) / svdParams{1}) + 1);% 2;
% end
if ~isfield(opts, 'scale') || ~isfield(opts.scale, varName)
    opts.scale.(varName) = [2 1];
end

model.restriction.(varName) = @restrictMatrixLeft;
model.prolongation.(varName) = @prolongateMatrixLeft;
model.operatorParams.(varName).prolongCoeff = 2;
model.operatorParams.(varName).normPColumns = false;
model.operatorParams.(varName).normPColumnsProlongationOnly = true;

Acoarse = restrictVar( model, A, opts, varName );

end

