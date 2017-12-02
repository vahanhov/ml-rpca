function varh = linearInterpolation( var, varLen, varh, varhLen, varargin )

%% Applies linear interpolation operator on a given variable to construct
%   its fine version. Used for prolongation.
% 
%   Author: Vahan Hovhannisyan, 2017.


if ~isempty(varargin)
    operParams = varargin{1};
    if isfield(operParams, 'normPColumns')
        normPColumns = operParams.normPColumns;
    else
        normPColumns = false;
    end
    %% THIS VIOLATES P = R^T, BUT PROVIDES WELL SCALED OPERATORS.
    if isfield(operParams, 'normPColumnsProlongationOnly')
        normPColumns = normPColumns || operParams.normPColumnsProlongationOnly;
    end
    if isfield(operParams, 'prolongCoeff')
        prolongCoeff = operParams.prolongCoeff;
    else
        prolongCoeff = 1;
    end
else
    normPColumns = false;
end

if normPColumns
    prev = var(:, 1);
else
    prev = 0;
end

for i = 1 : varLen
    curr = var(:, i); %getElementFromArray(var, i);
    next = curr;
    varhBlockPrev = 0.5 * curr;
    varh(:, 2 * i) = prolongCoeff * varhBlockPrev;
    %varh = setElementToArray(varh, 2 * i - 1, sigma * varhBlockPrev);
    varhBlockNext = 0.25 * (prev + next);
    varh(:, 2 * i - 1) = prolongCoeff * varhBlockNext;
    %varh = setElementToArray(varh, 2 * i, sigma * varhBlockNext);
    prev = next;
end
if varhLen > 2 * varLen
    varh(:, 2 * varLen + 1) = prolongCoeff * max(0.25, normPColumns * 0.5) * var(:, varLen);
else
    varh(:, 2 * varLen) = varh(:, 2 * varLen) + prolongCoeff * max(0, ~normPColumns * 0.25) * var(:, varLen);
end


end

