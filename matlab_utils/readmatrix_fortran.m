function matrix = readmatrix_fortran(filename)
    matrix = readmatrix(filename, delimitedTextImportOptions('DataLines', [1,Inf]));
    matrix = cellfun(@str2double, matrix);
    matrix = matrix(~isnan(matrix));
end