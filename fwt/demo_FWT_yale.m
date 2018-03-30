%% Run all facial shadow removal experiments
%
%   Author: Vahan Hovhannisyan, 2017.

data = 'CroppedYale/yaleB';

for i = 1:36
    if i == 14
        continue
    end
    demo_FWT([data num2str(i, '%02d')]);
end

