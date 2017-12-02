function [ R ] = getRestrictionMatrixRight( n, depth, varargin )
%% Returns the restirction matrix used
%   Used to test the restriction matrix, but not in the actual algorithms
%
%   Author: Vahan Hovhannisyan, 2017.

R = [];
for i = 1 : depth - 1
    if isempty(varargin)
        operator = @restrictMatrixRight;
        operParams = [];
    else
        operator = varargin{1};
        operParams = varargin{2};
    end

    I = eye(n);
    operParams.dimH = [n, floor(n / 2)];
    
    if isempty(R)
        R = operator(I, operParams);
    else
        R = R * operator(I, operParams);
    end
    n = floor(n / 2);
end

end

