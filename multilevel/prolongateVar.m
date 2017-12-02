function varsh = prolongateVar( var, model, opts, varName, varDims )
%% Applies the appropriate restriction operator on the variable up to the given depth
% 
%   Author: Vahan Hovhannisyan, 2016.

%% Iteratively prolongate vars to the finest level
varsh = var;
for level = opts.depth - 2 : -1 : 0
    prolongationOperator = model.prolongation.(varName);
    operParams = model.operatorParams.(varName);
    %operParams.imDims = floor(params.imDims ./ 2^level);
    % CAUTION! The following must be agreed with restriction
    operParams.dimh = max(floor(varDims ./...
        opts.scale.(varName) .^ level ), 1);
    varsh = prolongationOperator( varsh, operParams );
end

if isfield(model, 'postCoarse')
    varsh = model.postCoarse(coarserParams, varsh);
end

end
