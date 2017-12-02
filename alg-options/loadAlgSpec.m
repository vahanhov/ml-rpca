function [ struct ] = loadAlgSpec( struct, specName, varargin )
%% Automatically locates and loads relevant options and model files for 
%   each specific algorithm
%
%   Author: Vahan Hovhannisyan, 2017.

structSpec = str2func(specName);
structSpecFx = functions(structSpec);
if exist(structSpecFx.function)
    struct = structSpec(struct, varargin{:});
end

end

