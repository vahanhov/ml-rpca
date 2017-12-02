%% Runs all synthetic experiments
% 
%   Author: Vahan Hovhannisyan, 2017.


dims = {[50 100 100], [50 100 1000], [50 100 100], [50 100 1000]};
ranks = {2, 2, 5, 5};
maxtimes = {5, 10, 5, 10};

for i = 1:numel(dims)
    testRPCA('synth', dims{i}, ranks{i}, maxtimes{i});
end