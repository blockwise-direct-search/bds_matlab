function merge_pdf_order(outdir, outputfile, compdf_location)
% merge_pdf Merge pdf files in a folder.
% outdir: the folder where the pdf files are stored.
% outputfile: the full name of the merged pdf file (including ".pdf"). It
% shoud be the string of single quote.
% compdf_location: the location of the compdf file.

sort_order = {'no_noise', 'rotation', 'randomx0_10', 'rotation_noisy_guassian_-1',...
    'rotation_noisy_guassian_-2', 'rotation_noisy_guassian_-3', ...
    'rotation_noisy_guassian_-4', 'rotation_noisy_guassian_-5'};

cd(outdir);
% Initialize string variable.
pdfFiles = dir(fullfile(outdir, '*.pdf'));

% Store filename in a cell.
pdfNamesCell = cell(numel(pdfFiles), 1);
for i = 1:numel(pdfFiles)
    pdfNamesCell{i} = pdfFiles(i).name;
end

% Sort pdfNamesCell under the sort_order.
sortedIdx = zeros(length(sort_order), numel(pdfNamesCell));
for j = 1:numel(pdfNamesCell)
    for i = 1:length(sort_order)
        if contains(pdfNamesCell{j}, sort_order{i})
            sortedIdx(i, j) = 1;
        else
            sortedIdx(i, j) = 0;
        end
    end
end
[real_order, ~] = find(sortedIdx == 1);
pdfNamesCell = pdfNamesCell(real_order);

% Use the strjoin function to concatenate the elements in a cell array into a single string.
inputfiles = strjoin(pdfNamesCell, ' ');

% Remove spaces at the beginning of a string.
inputfiles = strtrim(inputfiles);

% Merge pdf.
system(['bash ', compdf_location, ' ', inputfiles, ' -o ', outputfile]);

end

