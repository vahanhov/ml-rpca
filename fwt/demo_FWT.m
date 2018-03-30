function demo_FWT(varargin)
%% Tests the FWT and ML-FWT methods on a given problem specified in varargin
%   Derived from the according function in the FW-T library
%
%   Author: Vahan Hovhannisyan, 2017.


currentPath = cd;

close all;

if isempty(varargin)
    debug = true;
    dbstop if error
    
	data = 'highway';
    % 'hall';
	% 'mall'
	% 'lobby'
	% 'hall'
    % 'copymachine'
else
    debug = false;
    
	data = varargin{1};
end

algorithms = {'FW_T', 'ML_FW_T'};
nAlgs = numel(algorithms);

delta = 1e-3;
rho = .75;  % sampling ratio; set to "1" for full observation

fprintf('**************************************************************\n')
fprintf(strcat(data, ' experiment', ' has started!'))

L0 = [];
S0 = [];
if ~strcmp(data, 'synth')
    imagePath = 'data/';
    
    if isdir(fullfile(imagePath, data))
        [ D, frameSize ] = loadImageFrames( imagePath, data, 5000 );
    elseif exist(strcat(imagePath, 'data_FWT/', data, '.mat'), 'file')
        load(strcat(imagePath, 'data_FWT/', data, '.mat'));
        %[ D, imDims, opts.L0, opts.S0, opts.GTidx ] = loadFile( userNameFull, maxImages, true );
    else
        warning('Data file does not exist!');
        return;
    end
    r = '';
else
    if numel(varargin) > 1 && ~isempty(varargin{2})
        r = varargin{2};
    else
        r = 7;
    end
    if numel(varargin) > 2 && ~isempty(varargin{3})
        sparsity = varargin{3};
    else
        sparsity = 0.5;
    end
    if numel(varargin) > 3 && ~isempty(varargin{4})
        incoh = varargin{4};
    else
        incoh = 1.25;%.25;
    end
    if numel(varargin) > 4 && ~isempty(varargin{5})
        dims = varargin{5};
    else
        dims = [20*r, 20*r, 100*r];
    end
    [ D, frameSize, L0, S0 ] = loadSynthetic( dims, r, incoh, sparsity );
end

par.GT.L = L0;
par.GT.S = S0;
clear L0 S0;

[m n] = size(D);

fprintf('data has been loaded: m = %d, n = %d; \n', m,n);

%% parameter tuning

if rho == 1
    
    fprintf('RPCA with full obseravation; \n');
    obs = D; Omega = ones(m,n);
    
else
    
    fprintf('RPCA with partial obseravation: ');
    Omega = rand(m,n)<=rho; % support of observation
    obs = Omega.*D; % measurements are made
    fprintf('observations are generated; \n');
    
end

% this is parameter to control noise level
% the smaller the noise, the smaller is delta
obs = obs/norm(obs, 'fro');
lambda_1 = delta*rho; 
lambda_2 = delta*sqrt(rho)/sqrt(max(m,n));

par.D = obs; 
par.lambda_1 = lambda_1; par.lambda_2 =lambda_2;
par.iter = 1000; 
par.epsilon = 1e-3; % stopping criterion
par.Omega = Omega;
par.showvideo = true; 
par.imDims = frameSize;
par.logs = debug;
par.maxtime = -1;

% display video or not
showvideo = debug;

% This is only for ML
% The less the noise, the deeper the levels can be
%% TODO: check the thoery to fix the deepest reasonable level
switch data
    case 'mall'
        par.opts.depth = 6;
    case 'hall'
        par.opts.depth = 5;
        %par.Daxtime = 2;
    case 'lobby'
        par.opts.depth = 6;
    case 'copymachine'
        par.opts.depth = 8;
        %par.epsilon = 1e-1;
    case 'synth'
        showvideo = false;
        depth = floor(log(n/7)/log(2)) + 1;
        par.opts.depth = depth;
        %par.epsilon = 1e-3;
    otherwise
        par.opts.depth = 4;
end

results = cell(nAlgs, 1);
times = cell(nAlgs, 1);

for algInd = 1 : nAlgs
    algName = algorithms{algInd};
    
    fprintf('**************************************************************\n')
    disp(['Let us try ', algName, ' method!']);
    fprintf('**************************************************************\n')
    
    alg = str2func(algName);
    par.starttime = cputime;
    output = alg(par); % main function
    time = cputime - par.starttime;
    
    decerr = norm(Omega .* (par.D - output.L - output.S), 'fro');
    output.nucnorm = norm(svd(output.L), 1);
    output.l1norm = norm(output.S, 1);
    output.objval = 0.5 * decerr^2 + par.lambda_1 * output.nucnorm ...
                    + par.lambda_2 * output.l1norm;
    output.decerr = decerr / norm(par.D, 'fro');
    
    results{algInd} = output;
    times{algInd} = time;
    
    % save results
    
    
    % summary of output
    fprintf('**************************************************************\n')
    fprintf('Summary of Output: \n')
    fprintf('%20s      %10s       %10s      %10s    %10s       %6s  %10s  \n',...
        'time', 'obj. val.', 'dec. err.', 'nuc. norm', 'l1 norm', 'rank', 'nnz');
    fprintf('%5s  %15.4d   %6d   %10d     %10d     %10d    %5d     %10d  \n', algName, time,...
        output.objval, output.decerr, output.nucnorm, output.l1norm,...
        rank(output.L), nnz(output.S)/numel(output.S));
    fprintf('**************************************************************\n')

    if showvideo && ~strcmp(data, 'synth')
        showVideo(par.D, frameSize, output);
    end
end

if ~debug
    file = num2str(now);
    if ~isempty(r)
        file = [file '_r' num2str(r) '_dim_n' num2str(dims(end))];
    end
    saveRegexp = '^(?!(par|D|GT|mask)$).';
    destRoot = fullfile(currentPath, 'fwt/results');
    save([destRoot '/exp_' data '_' file '.mat'], '-regexp', saveRegexp, '-v7.3');
    if ~strcmp(data, 'synth')
        %saveRPCA_samples( par, algorithms, results, data, destRoot, '', times, {}, 15);
    end
    close all;
end


if debug && isfield(results{1}.gt_error, 'L')
    col = {'b', 'r'};
    var = {'L', 'S'};
    for v = 1:numel(var)
        figure;
        for i = 1:nAlgs
            plot(results{i}.gt_error.(var{v}), col{i});
            hold on;
        end
        title(['Distance from ground truth: ', var{v}]);
        legend(algorithms);
        hold off;
    end
end

if debug
    error('Debug, you fool!');
end

end


