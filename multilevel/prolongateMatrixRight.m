function Ah = prolongateMatrixRight( A, operatorParams )
%% Applies n * n_h prolongation operator on a matrix of dimensions m * n.
% The prolongation operator is linear interpolation
%
%   Author: Vahan Hovhannisyan, 2017.

operParams.applyFromLeft = false;
Ah = prolongateMatrix( A, operatorParams );

end
