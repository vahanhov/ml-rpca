function [ D, imDims, L0, S0, GTidx ] = loadMatrix( data, varargin )
%% Loads experiments from a mat file
%   with ground truth data if present
%
%   Author: Vahan Hovhannisyan, 2017.

D = data.D;
imDims = data.imDims;

if isfield(data, 'GTidx')
    GTidx = data.GTidx;
else
    GTidx = ':';
end
if isfield(data, 'L0')
    L0 = data.L0;
else
    L0 = [];
end
if isfield(data, 'S0')
    S0 = data.S0;
else
    S0 = [];
end

if ~isempty(varargin)
    if varargin{1} < size(D, 2)
        nImages = varargin{1};
        D = D( :, 1 : nImages);
        L0 = L0( :, 1 : nImages);
        S0 = S0( :, 1 : nImages);
    end
end


end
