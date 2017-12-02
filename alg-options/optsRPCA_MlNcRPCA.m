function [ opts ] = optsRPCA_MlNcRPCA( opts, params, userName )
%% Sets parameters for ML-AltProj
%
%   Author: Vahan Hovhannisyan, 2017.


if ~isfield(opts, 'depth') || isempty(opts.depth)
    match = strsplit(userName, '/');
    match = match{1};
    switch match
        case 'highway'
            opts.depth = 2;
        case 'copymachine'
            opts.depth = 8;
        case 'walk'
            opts.depth = 6;
        case 'snow'
            opts.depth = 4;
        case 'gates'
            opts.depth = 8;
        case 'CroppedYale';
            opts.depth = 2;
            %opts.rank = 9;
        otherwise
            warning(['no optrions were specified for ' userName ', using defaults']);
            if strncmp(userName, 'CroppedYale/', numel('CroppedYale/'))
                opts.depth = 3;
            end
    end
end

opts = setDefaultOption(opts, 'depth', floor(log2(size(params.D, 2) / 25) + 1));
%opts.coarseThreshold = 1e-3;

end

