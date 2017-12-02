function [ y ] = soft( x, alpha )
%% Soft threasholding operator
%
%   Author: Vahan Hovhannisyan, 2016.

y = sign(x) .* max(abs(x) - alpha, 0);

end

