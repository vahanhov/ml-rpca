function varH = linearInterpolationTranspose( var, varLen, varH, varHLen, varargin )

%% Initial value of varH is given as input as it can be of different type,
%   e.g. matrix of zeros or an empty cell array.
% 
%   Author: Vahan Hovhannisyan, 2017.

if varHLen == varLen
    varH = var;
    return
end

if ~isempty(varargin)
    operParams = varargin{1};
    if isfield(operParams, 'rep')
        rep = operParams.rep;
    else
        rep = false;
    end
    if isfield(operParams, 'normPColumns')
        normPColumns = operParams.normPColumns;
    else
        normPColumns = false;
    end
    if isfield(operParams, 'restrictCoeff')
        restrictCoeff = operParams.restrictCoeff;
    else
        restrictCoeff = 1;
    end
else
    rep = false;
    normPColumns = false;
end

prev = var(:, 1);
if normPColumns
    prev = 2 * prev;
end
for i = 1 : varHLen % If scale.x > 2 concat as many restriction operators after each other as needed
    if rep
        ind = i : varHLen : varLen / 2; % Indices of starting points for each
    else
        ind = i;
    end
    %restriction operator to repeatedly use all elements of the fine
    %variable to construct the coarse variable, if rep == TRUE and scale > 2
    curr = 0; % Sum of every 2*ind-1 element in var
    next = 0; % Sum of every 2*ind element in var
    for j = 1 : length(ind)
        curr = curr + var(:, 2 * ind(j)); %getElementFromArray(var, 2 * ind(j) - 1);
        if 2 * ind(j) + 1 <= varLen
            next = next + var(:, 2 * ind(j) + 1); %getElementFromArray(var, 2 * ind(j));
        end
    end
    varHBlock =  0.5 * curr + 0.25 * (prev + next); % Weighted sum of three neighbouring columns
    varH(:, i) = restrictCoeff * varHBlock; % = setElementToArray(varH, i, varHBlock);
    prev = next;
end
if normPColumns
    varH(:, ind) = varH(:, ind) + 0.25 * next;
elseif next == 0
    varH(:, ind) = varH(:, ind) + 0.25 * curr;
end

end
