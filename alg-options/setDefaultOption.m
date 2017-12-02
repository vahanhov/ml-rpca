function [ opts ] = setDefaultOption( opts, option, value, varargin )
%% Sets a default value to option, if not already set.
%
%   Author: Vahan Hovhannisyan, 2017.

if ~isfield(opts, option)
    opts.(option) = value;
    
    if isnumeric(value) || islogical(value)
        valStr = num2str(value);
    else
        if isa(value, 'function_handle')
            valStr = func2str(value);
        else
            valStr = '';
        end
    end
    
    if isempty(varargin) || varargin{1}
        display(['No value is given for option ' option '. Setting to default ' valStr '.']);
    end
end

end

