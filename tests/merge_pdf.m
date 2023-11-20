function merge_pdf(outdir, outputfile, compdf_location)

cd(outdir);
% Initialize string variable.
pdfFiles = dir(fullfile(outdir, '*.pdf'));

% Store filename in a cell.
pdfNamesCell = cell(numel(pdfFiles), 1);
for i = 1:numel(pdfFiles)
    pdfNamesCell{i} = pdfFiles(i).name;
end

% Use the strjoin function to concatenate the elements in a cell array into a single string.
inputfiles = strjoin(pdfNamesCell, ' ');

% Remove spaces at the beginning of a string.
inputfiles = strtrim(inputfiles);

% Merge pdf.
% outputfile = 'merged_CBDS_5_200_big_plain.pdf';
%system(['bash ', '/home/lhtian97/Documents/bds/tests/private/compdf', ' ', inputfiles, ' -o ', outputfile]);
keyboard
system(['bash ', compdf_location, ' ', inputfiles, ' -o ', outputfile]);


end

