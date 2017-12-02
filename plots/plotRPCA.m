function plotRPCA( imDims, varargin )
%% Creates plots from the experimental results
%
%   Author: Vahan Hovhannisyan, 2017.

if ~exist('imshow')
    warning('imshow was not found!');
    return
end

layout.xI = 20;
layout.yI = 20;
layout.gap = 2;
layout.gap2 = 1;

% Hardcode names!
names = {'Original', 'Low Rank', 'Sparse'};

for i = 1 : numel(varargin)
    plotResults(varargin{i}, imDims, names{i}, layout);
end

end




function plotResults(dataMatrix, imDims, plotTitle, layout)

if nargin < 4
    xI = ceil(sqrt(numImage));
    yI = ceil(numImage/xI);

    gap = 2;
    gap2 = 1;
else
    xI = layout.xI;
    yI = layout.yI;

    gap = layout.gap;
    gap2 = layout.gap2;
end

container = ones(imDims(1) + gap, imDims(2) + gap);
bigpic = cell(xI,yI);

[m, n] = size(dataMatrix);

for i = 1 : xI
    for j = 1 : yI
        if yI * (i - 1) + j > n
            bigpic{i,j} = ones(imDims(1) + gap, imDims(2) + gap);
        else
            container((gap2 + 1) : (end - gap2), (gap2 + 1) : (end - gap2))...
                            = reshape(dataMatrix(:, yI * (i - 1) + j), imDims);
            bigpic{i,j} = container;
        end
    end
end
figure
imshow(cell2mat(bigpic), [], 'Border', 'tight')
title(plotTitle);


end





