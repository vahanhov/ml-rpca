function Ah = prolongateMatrix( A, operParams )
%% Applies n * n_H prolongation operator on a m * n matrix from left,
%   or its transpose from right.
% 
%   Author: Vahan Hovhannisyan, 2017.

if isfield(operParams, 'applyFromLeft')
    applyFromLeft = operParams.applyFromLeft;
else
    applyFromLeft = false;
end

if isfield(operParams, 'prolongation')
    operator = operParams.prolongation;
else
    operator = @linearInterpolation;
end

dimh = operParams.dimh;

if applyFromLeft
    A = A';
    dimh = fliplr(dimh);
end

[m, n] = size(A);
if all(dimh == [m n])
    Ah = A;
else % Apply interpolation type restriction operator taking each column as a vector entry
    Ah = zeros(dimh);
    Ah = operator( A, n, Ah, dimh(2), operParams );
end

if applyFromLeft
    Ah = Ah';
end

% normalize columns.
if isfield(operParams, 'prolConst')
    Ah = Ah .* operParams.prolConst;
end


end

