%% Runs all real experiments
%   Change test to 'face' or 'video' for the according experiment
% 
%   Author: Vahan Hovhannisyan, 2017.

clear
close all

test = 'face'; %'face' 'video'

if strcmp(test, 'face')
    [dicts{1:39}] = deal('CroppedYale');
    for i = 1 : 9
        dicts{i} = [dicts{i} '/yaleB0' num2str(i)];
    end
    for i = 10 : 39
        dicts{i} = [dicts{i} '/yaleB' num2str(i)];
    end
    testind = [1, 2, 10];
    for i = 1 : numel(dicts(testind))
        testRPCA(dicts{testind(i)});
        close all;
    end
elseif strcmp(test, 'video')
    dicts = {'highway', 'copymachine', 'walk', 'gates'};
    
    for i = 1 : numel(dicts)
        testRPCA(dicts{i});
        close all;
    end
    
end



