function Ah = prolongateMatrixLeft(A, operParams)
%% Applies n * n_h prolongation operator on a matrix of dimensions m * n.
% The prolongation operator is linear interpolation
%
%   Author: Vahan Hovhannisyan, 2017.

operParams.applyFromLeft = true;
Ah = prolongateMatrix( A, operParams );


end

