function [ D, imDims, L0, S0, GTidx ] = loadTensor( data, varargin )
%% Loads data from a mat file stored a tensor of w * h * n dimensions
%
%   Author: Vahan Hovhannisyan, 2017.


dataX = data.X;

if isfield(data, 'GTidx')
    GTidx = data.GTidx;
else
    GTidx = ':';
end
if isfield(data, 'GT')
    dataGT = data.GT;
end
if isfield(data, 'L0')
    dataL0 = data.L0;
end
if isfield(data, 'S0')
    dataS0 = data.S0;
end

if ~isempty(varargin)
    if varargin{1} < size(dataX, 3)
        nImages = varargin{1};
        dataX = dataX( :, :, 1 : nImages);
        dataGT = dataGT( :, :, 1 : nImages);
        dataL0 = dataL0( :, :, 1 : nImages);
        dataS0 = dataS0( :, :, 1 : nImages);
    else
        nImages = size(dataX, 3);
    end
end

imDims = size(dataX(:, :, 1));
D = reshape(dataX, [prod(imDims), nImages]);
if isfield(data, 'L0')
    L0 = reshape(dataL0, [prod(imDims), size(dataL0, 3)]);
else
    L0 = [];
end
if isfield(data, 'S0')
    S0 = reshape(dataS0, [prod(imDims), size(dataS0, 3)]);
else
    S0 = [];
end


end

