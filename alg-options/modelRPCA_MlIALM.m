function [ model ] = modelRPCA_MlIALM( model, params )
%% Sets methods and their parameters for ML-IALM
%
%   Author: Vahan Hovhannisyan, 2017.


model.restriction.L = @restrictMatrixRight;
model.operatorParams.L.prolongCoeff = 2;
model.operatorParams.L.restrictCoeff = 1;
model.operatorParams.L.normPColumns = false;
model.operatorParams.L.normPColumnsProlongationOnly = true;

end

