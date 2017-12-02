function [ opts ] = optsRPCA_MlIALM( opts, params, dictionary )
%% Sets up parameters for the MLIALM algorithm for each experiments accordingly
%
%   Author: Vahan Hovhannisyan, 2017.

opts = optsRPCA_IALM( opts, params, dictionary );
match = strsplit(dictionary, '/');
match = match{1};
switch match
    case 'highway'
        opts.depth = 2;
    case 'copymachine'
        opts.depth = 4;
    case 'walk'
        opts.depth = 4;
        %opts.rho = 1.5;
        %opts.epsRound = 1e-5;
    case 'gates'
        opts.depth = 4;
        opts.rho = 1.3;
        opts.muCoeff = 20;
    case 'synth'
        opts.depth = 3; %floor(log(size(params.D, 2)/7)/log(2));% 3;
    case 'CroppedYale';
        opts.depth = 3;
        opts.muCoeff = 50;
        opts.rho = 1.2;
        opts.lambdaMlCoeff = 5;
        opts.epsRound = 1e-5;
        if strcmp(dictionary, 'CroppedYale/yaleB02')
            opts.lambdaMlCoeff = 8;
        end
    otherwise
        warning(['no optrions were specified for ' userName ', using defaults']);
        if strncmp(dictionary, 'CroppedYale/', numel('CroppedYale/'))
            opts.depth = 3;
        end
end

opts = setDefaultOption(opts, 'depth', floor(log2(size(params.D, 2) / 25) - 1));
if opts.depth <= 1
    error(['Doesnt make sense to use depth ' num2str(opts.depth)]);
end
opts.coarseVarNames = {'L'};

% R = getRestrictionMatrixRight(size(params.D, 2), opts.depth);
opts.svd_coarse_coeff = 1; %norm(R, 'inf'); % <1 values make the thresholding more forgiving

end
