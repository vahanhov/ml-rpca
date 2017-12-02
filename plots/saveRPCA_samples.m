function saveRPCA_samples( params, methods, xEsts, dic, destRoot, r,...
    times, opts, varargin )
%% Saves video frame samples from experimental results
%
%   Author: Vahan Hovhannisyan, 2017.

%% Set up
destDir = fullfile(destRoot, 'img' );
if ~exist(destDir, 'dir')
    mkdir(destRoot, 'img');
end
ext = 'eps';
if isempty(varargin)
    switch dic
        case 'gates'
            imInd = 220; %round(size(params.D, 2) / 2);
        case 'CroppedYale/yaleB02'
            imInd = 13;
        otherwise
            imInd = round(size(params.D, 2) / 2);
    end
else
    imInd = varargin{1};
end

iptsetpref('ImshowBorder','tight');

%% Original image
orig = reshape(params.D(:, imInd), params.imDims);
h = figure;
imshow(orig, []);
set(gca,'position',[0 0 1 1],'units','normalized');

print(h, ['-d', ext, 'c'], [destDir '/' dic '_orig' r]);

for i = 1 : numel(xEsts)
    if isempty(xEsts{i})
        continue
    end
    %% Low rank sample
    lowrank = reshape(xEsts{i}.L(:, imInd), params.imDims);
    h = figure;
    imshow( lowrank, [] );
    set(gca,'position',[0 0 1 1],'units','normalized');
    print(h, ['-d', ext, 'c'], [destDir '/' dic '_lowrank_' methods{i} r]);
    %% Sparse sample
    sparse = reshape(xEsts{i}.S(:, imInd), params.imDims);
    h = figure;
    imshow( sparse, [] );
    set(gca,'position',[0 0 1 1],'units','normalized');
    print(h, ['-d', ext, 'c'], [destDir '/' dic '_sparse_' methods{i} r]);
    %% Further crucial data
    ranks{i} = rank(xEsts{i}.L);
    nucNorm{i} = norm(svd(xEsts{i}.L, 'econ'), 1);
    l1Norm{i} = norm(xEsts{i}.S, 1);
    decError{i} = num2str(norm(params.D - xEsts{i}.L - xEsts{i}.S, 'fro'));
end

if numel(methods) == 2
    alg = ['_' methods{2}];
else
    alg = '';
end
for i = 1 : numel(opts)
    if isfield(opts{i}, 'L0')
        opts{i} = rmfield(opts{i}, 'L0');
    end
    if isfield(opts{i}, 'S0')
        opts{i} = rmfield(opts{i}, 'S0');
    end
end
save([destDir '/' dic '_data' alg r '.mat'], '-v7.3', 'ranks', 'nucNorm', 'l1Norm', 'decError', 'methods', 'dic', 'times', 'opts');

end

