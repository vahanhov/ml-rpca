function [ Aapprox, U, ev, V ] = mlLinNuc( A, Omega, opts )
%% Calculates a rank one approximation for the given matrix A 
%       using using its coarse model
% A                 given matrix to be approximated
% Omega             the projection operator, if empty then identity matrix
% opts              options for restriction and prolongation 
%
%   Author: Vahan Hovhannisyan, 2017.

varName = 'A';
[ ACoarse, model, opts ] = restrictVarRPCA_left( A, varName, opts );

coeff = sqrt(size(A, 1) / size(ACoarse, 1));

if ~isempty(Omega)
    ACoarse = Omega .* ACoarse;
end

[ U, s, V ] = lansvd( ACoarse, 1, 'L' );
ev = coeff .* s;

ACoarseApprox = -U * V'; % Scale according to the restriction operator??
varDims = size(A); %% TODO: prolongate only V
Aapprox = prolongateVar( ACoarseApprox, model, opts, varName, varDims ) ./ coeff;

% params.inits.Vt = zeros(1, size(params.inits.(varName), 2));
% model.prolongation.Vt = model.prolongation.(varName);
% model.operatorParams.Vt = model.operatorParams.(varName);
% opts.scale.Vt = opts.scale.(varName);
% Vht = prolongateVar( model, V', params, opts, 'Vt' );
% Aapprox = -U * Vht;

end

