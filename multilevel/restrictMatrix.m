function AH = restrictMatrix( A, operParams )
%% Applies n * n_H restriction operator on a m * n matrix from left,
%   or its transpose from right.
% 
%   Author: Vahan Hovhannisyan, 2016.

if isfield(operParams, 'applyFromLeft')
    applyFromLeft = operParams.applyFromLeft;
else
    applyFromLeft = false;
end

if isfield(operParams, 'restriction')
    operator = operParams.restriction;
else
    operator = @linearInterpolationTranspose;
end


dimH = operParams.dimH;

if applyFromLeft
    A = A';
    dimH = fliplr(dimH);
end

[m, n] = size(A);
if n > 1
    if all(dimH == [m n])
        AH = A;
    else % Apply interpolation type restriction operator taking each column as a vector entry
        zerosColumn = zeros(m, 1);
        AH = operator( A, n, zerosColumn, dimH(2), operParams );
    end
else
    AH = A;
end

if applyFromLeft
    AH = AH';
end

% normalize rows.
if isfield(operParams, 'restConst')
    AH = AH .* operParams.restConst;
end

end

