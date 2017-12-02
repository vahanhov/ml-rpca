%% Script to test restriction and prolongation operators
%
%   Use this for new operators to make sure they are transpose to each
%   other
%
%   Author: Vahan Hovhannisyan, 2017.

clear

imDims = [5, 3];
n = prod(imDims);

nH = floor( n / 2 );
Ih = eye(n, n);
IH = eye(nH, nH);

operParams.imDims = imDims;
operParams.dimH = [n, nH];
operParams.dimh = [nH, n];
operParams.prolongCoeff = 1;
operParams.normPColumns = true;


R = restrictMatrixRight( Ih, operParams );
RT = R';
P = prolongateMatrixRight( IH, operParams );

if ~all(all(operParams.prolongCoeff * RT == P))
    warning('-- ERROR --!');
    disp(num2str(RT - P));
else
    disp('-- CONGRATULATIONS! --');
end