function merge_pdf(outdir, outputfile, compdf_location)
% merge_pdf Merge pdf files in a folder.
% outdir: the folder where the pdf files are stored.
% outputfile: the full name of the merged pdf file (including ".pdf").
% compdf_location: the location of the compdf file.

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
system(['bash ', compdf_location, ' ', inputfiles, ' -o ', outputfile]);

end

