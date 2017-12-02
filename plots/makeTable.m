function makeTable( indir, outfile )
%% Creates and saves a latex table from a .mat file containing experimental results
%
%   Author: Vahan Hovhannisyan, 2017.

files = list_files(indir);

tabledata = {''};

for i = 1:numel(files)
    f = files{i};
    [~, fname, fext] = fileparts(f);
    if fext ~= 'm'
        continue
    end
    
    load(f);
    
    tabledata{i, 1} = strtok(fname, '_');
    if ~exist('algorithms', 'var')
        continue
    end
    nalgs = numel(algorithms);
    for j = 1:nalgs
        tabledata{i, 5*(j-1)+(nalgs-j+1)} = num2str(times{j}, 2);
        tabledata{i, 5*(j-1)+(nalgs-j+2)} = num2str(results{j}.objval, 2);
        tabledata{i, 5*(j-1)+(nalgs-j+3)} = num2str(rank(results{j}.L), 2);
        tabledata{i, 5*(j-1)+(nalgs-j+4)} = num2str(nnz(results{j}.S) / (m*n), 2);
    end
end

inputlatex.data = tabledata;
latexdata = char(latexTable(inputlatex));
save(fullfile(indir, outfile), 'latexdata');

end

