function run_pdfjam(outdir, num_feature, pdfname)

cd(outdir);

% Check if the file name ends with .pdf.
if ~endsWith(pdfname, '.pdf')
    % If not, add the extension.
    pdfname_output = [pdfname, '.pdf'];
else
    % If yes, use the same name.
    pdfname_output = pdfname;
end

epsname_output = strrep(pdfname_output, '.pdf', '.eps');

% Create pdfjam command
pdfjamCommand = sprintf('pdfjam "%s" --nup 1x%d --fitpaper true --outfile "%s"', pdfname, num_feature, pdfname_output);

% Execute the system command
[status, cmdout] = system(pdfjamCommand);

% Check the status of the command
if status == 0
    fprintf('Output file generated successfully: %s\n', pdfname_output);
else
    fprintf('Command execution failed:\n%s\n', cmdout);
end

% Execute pdfcrop command to cut the white margins
pdfcropCommand = sprintf('pdfcrop "%s" "%s"', pdfname_output, pdfname_output);
[status1, cmdout1] = system(pdfcropCommand);

% Check the status of the command
if status1 == 0
    fprintf('Crop PDF successfully: %s\n', pdfname_output);
else
    fprintf('pdfcrop command execution failed:\n%s\n', cmdout1);
end

% Convert the cropped PDF to EPS format
pdftopsCommand = sprintf('pdftops -eps "%s" "%s"', pdfname_output, epsname_output);
[status2, cmdout2] = system(pdftopsCommand);

% Check the status of cropping the PDF to EPS format
if status2 == 0
    fprintf('Successfully generated EPS file: %s\n', epsname_output);
else
    fprintf('pdftops command execution failed:\n%s\n', cmdout2);
end

end

