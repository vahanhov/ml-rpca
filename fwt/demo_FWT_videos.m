%% Run all video background extraction experiments
%
%   Author: Vahan Hovhannisyan, 2017.

clear

data = {'highway', 'hall', 'mall', 'lobby', 'copymachine'};

for i = 1 : numel(data)
    demo_FWT(data{i});
end
