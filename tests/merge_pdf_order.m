function merge_pdf_order(outdir, outputfile, compdf_location)
% merge_pdf Merge pdf files in a folder.
% outdir: the folder where the pdf files are stored.
% outputfile: the full name of the merged pdf file (including ".pdf"). It
% should be the string of single quote.
% compdf_location: the location of the compdf file.

sort_order = {'plain', 'rotation_no_noise', 'randomx0_10', 'rotation_noisy_gaussian_-1',...
    'rotation_noisy_gaussian_-2', 'rotation_noisy_gaussian_-3', ...
    'rotation_noisy_gaussian_-4', 'rotation_noisy_gaussian_-5', ...
    'rotation_noisy_gaussian_-6', 'rotation_noisy_gaussian_-7', ...
    'rotation_noisy_gaussian_-8', "-1_no_rotation", "-2_no_rotation", ...
    "-3_no_rotation", "-4_no_rotation", "-5_no_rotation", ...
    "-6_no_rotation", "-7_no_rotation", "-8_no_rotation"};

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
% Since the element in the real_order may be greater than the maximum index of sortedIdx,
% We need to use the index of the real_order according to the sortedIdx to
% sort the pdfNamesCell.
[~, real_order_index] = sort(real_order);
keyboard
pdfNamesCell = pdfNamesCell(real_order_index);
% Use the strjoin function to concatenate the elements in a cell array into a single string.
inputfiles = strjoin(pdfNamesCell, ' ');

% Remove spaces at the beginning of a string.
inputfiles = strtrim(inputfiles);

% Merge pdf.
system(['bash ', compdf_location, ' ', inputfiles, ' -o ', outputfile]);

end

