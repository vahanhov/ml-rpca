function varsH = restrictVar( model, var, opts, varName )
%% Applies the appropriate restriction operator the variable up to the specified depth
% 
%   Author: Vahan Hovhannisyan, 2017.

%% Iteratively restrict vars to the coarsest level
varsH = var;
for level = 1 : opts.depth - 1
    restrictionOperator = model.restriction.(varName);
    operParams = model.operatorParams.(varName);
    %operParams.imDims = floor(params.imDims ./ 2^(level-1));
    % CAUTION! The following must be agreed with prolongation
    operParams.dimH = max(floor(size(var) ./...
        opts.scale.(varName) .^ level ), 1);
    varsH = restrictionOperator( varsH, operParams );
end

if any(size(varsH) <= 1 & size(var) == 1)
    error('Too many levels! Too many levels!');
end

end

