function testRPCA(varargin)
%% Tests RPCA methods on various datasets
% 
%   Author: Vahan Hovhannisyan, 2017.

dbstop if error
close all

%% Choose the methods to test
%methodNames = {'NcRPCA', 'MlNcRPCA'};
methodNames = {'IALM', 'MlIALM'};


if ~isempty(varargin)
    debug = false;
else
    debug = true;
end

currentPath = cd;

% input path
if isempty(varargin)
    userName = 'highway'; %'CroppedYale/yaleB01'; % 'copymachine'; 'highway'; 'makeup'; 'snow'; 'gates'; 'walk'; 'yaleB01'; 'yaleB11'; 'YaleB96x84'; 'synth'
else
    userName = varargin{1};
end

%% Experiment options
opts.lambdaCoeff = 1;
maxImages = 5000; % Maximum number of images to load.
% maxTime is only for IALM methods
match = strsplit(userName, '/');
match = match{1};
switch match
    case 'highway'
        opts.maxTime = 2;
        params.rank = 1;
        opts.lambdaCoeff = 0.5;
        opts.muCoeff = 5;
    case 'copymachine'
        opts.maxTime = 50;
        params.rank = 1;
        opts.lambdaCoeff = 0.2;
        opts.muCoeff = 5;
    case 'walk'
        opts.maxTime = 50;
        params.rank = 1;
        opts.lambdaCoeff = 0.4;
        opts.muCoeff = 10;
    case 'gates'
        opts.maxTime = 1000;
        params.rank = 2;
        opts.lambdaCoeff = 1;
        opts.muCoeff = 10;
    case 'synth'
        dims = varargin{2};
        params.rank = varargin{3};
        opts.maxTime = varargin{4};
    case 'FR' % Only for debugging
        opts.maxTime = 100;
        params.rank = 20;
        maxImages = 100;
    case 'CroppedYale'
        opts.maxTime = 5;
        params.rank = 9;
        if strcmp(userName, 'CroppedYale/yaleB02')
            opts.lambdaCoeff = 0.5;
            opts.incoh = 1;
        end
    otherwise % Used for eavery CroppedYale database
        warning(['no optrions were specified for ' userName ', using defaults']);
        opts.maxTime = 100;
        params.rank = 10;
end

robust = true;
params.robust = robust;

if ispc
    imagePath = fullfile(currentPath, 'data_win');
else
    imagePath = fullfile(currentPath, 'data');
end

if maxImages > 0
    disp(['Using at most the first ' num2str(maxImages) ' images.']);
end

%% Get training images
if strcmp( userName, 'FR' )
    dicPath = '\\fs-vol-bitbucket.doc.ic.ac.uk\bitbucket\vh13\FR_exp\';
    [ D, ~, ~, ~, imDims ] = loadDictionaryFR( 1, dicPath );
elseif strcmp( userName, 'synth' )
    [ D, imDims, opts.L0, opts.S0 ] = loadSynthetic( dims, params.rank );
else
    userNameFull = fullfile(imagePath, userName);
    if isdir(userNameFull)
        [ D, imDims ] = loadImageFrames( imagePath, userName, maxImages );
    else
        [ D, imDims, opts.L0, opts.S0, opts.GTidx ] = loadFile( userNameFull, maxImages, true );
        %D = D ./ norm(D);
    end
end
opts = setDefaultOption(opts, 'GTidx', ':');

saveRegexp = '^(?!(params|D|GT|mask)$).';

%% Set problem params
params.imDims = imDims;
params.D = D;

[m, n] = size(D);
initCoef = 1;%1e-3;
rng(10);
params.inits.L = initCoef * rand(m, n); %zeros(m, n);
if robust
    params.inits.S = initCoef * rand(m, n); %zeros(m, n);
end

params.nu = opts.lambdaCoeff / sqrt(max(size(D))); %% Required by theory. Do not change.
opts.maxIter = 1e5;
opts.haltTime = 1e4;
opts.tolerance = 1e-7;

opts.logs = true;
if maxImages >= 500
    opts.logInterval = 10;
else
    opts.logInterval = 10;
end

opts.timeStepLen = opts.maxTime / 20; % measure time at 20 spots

%% Run the algorithm
file = num2str(now);

numMethods = numel(methodNames);
results = cell(numMethods, 1);
fEsts = cell(numMethods, 1);
nIters = cell(numMethods, 1);
logData = cell(numMethods, 1);
times = cell(numMethods, 1);
output = cell(numMethods, 1);

% Store algorithm specific options for logging purposes
optsAlgs = cell(numMethods, 1);
modelAlgs = cell(numMethods, 1);

for methodIndex = 1 : numMethods
    fprintf(['\n\n Running with ' methodNames{methodIndex}...
                    ' on a ' num2str(m) ' x ' num2str(n) ' database.\n']);
    algorithm = str2func(['Solve' methodNames{methodIndex}]);
    
    optsAlgs{methodIndex} = loadAlgSpec( opts,...
        ['optsRPCA_' methodNames{methodIndex}], params, userName );
    
    t0 = cputime;
    opts.startTime = t0;
    [ xStar, FStar, nIter, logDataAlg ] = algorithm( modelAlgs{methodIndex},...
        params, optsAlgs{methodIndex} );
	tEst = cputime - t0;

    %% show results
    disp(['CPU time: ' num2str(tEst) '; nIter = ' num2str(nIter)]);
    
    %% Log results
    results{methodIndex} = xStar;
    nIters{methodIndex} = nIter;
    times{methodIndex} = tEst;
    logData{methodIndex} = logDataAlg;
    
    %% Print results
    if isfield(opts, 'S0') && ~isempty(opts.S0)
        ref = opts.S0;
    else
        ref = params.D;
    end
    
    if robust
        %mask = soft(xStar.L - params.D, 1e-6);
        if debug && ~strcmp(userName, 'synth')
            plotRPCA( params.imDims, ref, xStar.L, xStar.S );%, xEst.augLagrVar );
        end
        l1Norm = norm(xStar.S, 1);
        sparsity = nnz(xStar.S) / numel(params.D);
        decError = norm(params.D - xStar.L - xStar.S, 'fro');
    else
        if ~strcmp(userName, 'synth')
            plotRPCA( params.imDims, params.D, xStar.L );
        end
        l1Norm = nan;
        sparsity = nan;
        decError = norm(params.D - xStar.L, 'fro')
    end
    decError = decError / norm(D, 'fro');
    nucNorm = norm(svd(xStar.L, 'econ'), 1);
    rankL = rank(xStar.L);
    
    disp(['Rank( L ) = ' num2str(rankL)...
        '; || L ||_* = ' num2str(nucNorm)...
        '; || S ||_0 = ' num2str(sparsity)...
        '; || S ||_1 = ' num2str(l1Norm)...
        '; dec error = ' num2str(decError)]);
    
    if isempty(FStar) || isnan(FStar)
        FStar = nucNorm + params.nu * l1Norm;
    end
    fEsts{methodIndex} = FStar;
    output{methodIndex}.rank = rankL;
    output{methodIndex}.nucNorm = nucNorm;
    output{methodIndex}.sparsity = sparsity;
    output{methodIndex}.l1Norm = l1Norm;
    output{methodIndex}.decError = decError;
    
    if ~debug
        if isfield(params, 'rank') && strcmp(userName, 'synth')
            r = ['_r' num2str(params.rank)];
        else
            r = '';
        end
        if ~ispc
            save([fullfile(currentPath, 'rpca/results'), '/exp_' userName '_' file r '.mat'], '-regexp', saveRegexp, '-v7.3');
        end
    else
        r = '';
    end
end

if ~ispc
    % output
    destRoot = fullfile(currentPath, 'rpca/results');
    destDir = fullfile(destRoot, 'plots' );
    if ~exist(destDir, 'dir')
        mkdir(destRoot, 'plots');
    end
    if (isfield(opts, 'L0') && ~isempty(opts.L0)) || (isfield(opts, 'S0') && ~isempty(opts.S0))
        plotRPCA_GT( params, methodNames, results, logData, [userName '_' file r], destRoot );
    end
    if ~debug && ~strcmp(userName, 'synth')
        saveRPCA_samples( params, methodNames, results, userName, destRoot, r, times, optsAlgs, 15);
        close all
    end
end
if debug
    error('waiting for debug analysis');
end


end