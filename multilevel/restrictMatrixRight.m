function AH = restrictMatrixRight( A, operParams )
%% Applies n * n_h prolongation operator on a matrix of dimensions m * n.
% The prolongation operator is linear interpolation
% 
%   Author: Vahan Hovhannisyan, 2016.

operParams.applyFromLeft = false;
AH = restrictMatrix( A, operParams );

end

