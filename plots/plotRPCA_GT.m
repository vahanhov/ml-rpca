function plotRPCA_GT( params, algs, results, logData, dic, destRoot )
%% Creates plots from experimental results containing ground truth data
%
%   Author: Vahan Hovhannisyan, 2017.

destDir = fullfile(destRoot, 'plots' );
if ~exist(destDir, 'dir')
    mkdir(destRoot, 'plots');
end
ext = 'eps';

attrNames = {'errL', 'errS', 'frob_err'};
measures = {'relative error from ground truth',...
    'RE from ground truth',...
    'RRE'};
rangeIters = {0, 0, 1};
logVals = {0, 0, 1};

for attrInd = 1 : numel(attrNames)
    h = plotGTAttr( logData, attrNames{attrInd}, algs, measures{attrInd},...
        rangeIters{attrInd}, logVals{attrInd} );
    if ~isempty(h)
        pause(1);
        print(h, ['-d', ext, 'c'], [destDir '/' dic '_' attrNames{attrInd} '']);
    end
end

end



function h = plotGTAttr( logData, attrName, algs, measure, rangeIters, logVal )

allmethods = {'IALM', 'MlIALM', 'NcRPCA', 'MlNcRPCA'};
plotColors = containers.Map(allmethods, {'r', 'k', 'b', 'm'});
markers = containers.Map(allmethods, {'o', '+', '*', '.'});
plotStyles = containers.Map(allmethods, {'-', '--', ':', '-.'});
algNames = containers.Map(allmethods, {'IALM', 'ML-IALM', 'AltProj', 'ML-AltProj'});
fs = 28;

h = [];
for algInd = 1 : numel(logData)
    logDataAlg = logData{algInd};
    if ~isfield(logDataAlg, attrName)
        warning('The requested attribute had not been stored.');
        break;
    end
    if algInd == 1
        h = figure;
    end
    
    attr = logDataAlg.(attrName);
    if logVal
        attr = log(attr);
        measureAlg = [measure ' (log)'];
    else
        measureAlg = measure;
    end
    
    if rangeIters
        range = 1 : numel(attr);
        rangeLabel = 'iterations';
    else
        range = logDataAlg.timeSteps;
        rangeLabel = 'CPU time (seconds)';
    end
    
    plot(range, attr,...
        'Color', plotColors(algs{algInd}),...
        'Marker', markers(algs{algInd}),...
        'MarkerSize', 18);
    hold on    
end

if ~isempty(h)
    set(gca, 'FontSize', fs);
    l = legend(values(algNames, algs));
    set(l, 'Interpreter', 'latex', 'fontsize', fs);
    xlabel(rangeLabel, 'Interpreter', 'latex', 'fontsize', fs);
    ylabel(measureAlg, 'Interpreter', 'latex', 'fontsize', fs);
    %title(attrName, 'fontsize', fs);
end

end


