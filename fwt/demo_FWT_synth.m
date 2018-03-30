%% Run all synth experiments
%
%   Author: Vahan Hovhannisyan, 2017.


clear
dims = {[100 100 500], [100 100 5000], [100 100 500], [100 100 5000]};
ranks = {2, 2, 5, 5};

for i = 1 : numel(ranks)
    demo_FWT('synth', ranks{i}, '', '', dims{i});
end
