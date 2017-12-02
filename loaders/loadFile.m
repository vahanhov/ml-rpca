function [ D, imDims, L0, S0, GTidx ] = loadFile( fileName, varargin )
%%  Loads an experiment file from a file
%   with ground truth data if present
%
%   Author: Vahan Hovhannisyan, 2017.

data = load(fileName);

if isfield(data, 'X')
    [ D, imDims, L0, S0, GTidx ] = loadTensor(data, varargin{:});
elseif isfield(data, 'D')
    [ D, imDims, L0, S0, GTidx ] = loadMatrix(data, varargin{:});
end


% Normalize if required
if numel(varargin) > 1
    if varargin{2}
        for i = 1 : size(D, 2)
            D(:, i) = D(:, i) ./ max(D(:, i)); %norm(imcol);
%             if ~isempty(GT)
%                 GT(:, i) = GT(:, i) ./ max(GT(:, i)); %norm(imcol);
%             end
        end
    end
end


disp(['Loaded ' num2str(size(D, 2)) ' images from ' num2str(fileName) '.']);

end

