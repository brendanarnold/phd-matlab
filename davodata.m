delimiter = MiscFns.tab;
header_line = 1;
input_extension = '.*_(.*)-.*_spliced.dat$';
output_extension = '.avg';
x_col = 2;
x_col_header = 'B(T)';
y_col = 3;
bin_size = 0.2;

% Say some stuff
disp('davodata.m');
disp('  Averages data over the specified bin sizes (set in file)');
disp('  similar to the Pascal program of the same name.');
disp(sprintf('Bin size: %.2f', bin_size));
disp(sprintf('X column: %d', x_col));
disp(sprintf('Y column: %d', y_col));

file_info = dir();
for index = 1:length(file_info)
    filename = file_info(index).name; % Cache the filename
    if regexp(filename, input_extension)
        % Compile the column headers
        [y_col_header] = regexp(filename, input_extension, 'tokens');
        % Why are the tokens stored within a cell within a cell?!
        y_col_header = y_col_header{1};
        y_col_header = y_col_header{1};
        y_col_header = regexprep(y_col_header, 'p', '.');
        colheaders = {x_col_header y_col_header 'StdevY'};
        file = MiscFns.import_csv(file_info(index).name, delimiter, header_line);
        
        % DO SOMETHING TO EACH FILE HERE
        
        % Replicates the 'Davodata' program in MATLAB, tested output matches
        % exactly.
               
        max_x = max(file.data(:, x_col));
        min_x = min(file.data(:, x_col));
        [num_points, num_cols] = size(file.data);
        num_bins = floor((max_x - min_x)/bin_size) + 1;
        % Initialise some matrices for performance
        bins = zeros(num_bins, 4);
        output = zeros(num_bins, 3);
        
        % There is a probably a more 'MATLAB' way to do this, but this is not
        % a big deal..
        
        % Bin the data
        for line = file.data'
            x = line(x_col); 
            y = line(y_col);
            bin_index = floor((x-min_x)/bin_size) + 1;
            bins(bin_index, 1) = bins(bin_index, 1) + x;
            bins(bin_index, 2) = bins(bin_index, 2) + y;
            bins(bin_index, 3) = bins(bin_index, 3) + y*y;
            bins(bin_index, 4) = bins(bin_index, 4) + 1;
        end
        
        % Calculate the average x, y and stdev of y
        counter = 1;
        for bin = bins'
            total_x = bin(1);
            total_y = bin(2);
            total_y2 = bin(3);
            tally = bin(4);
            if tally > 0
                avg_x = total_x/tally;
                avg_y = total_y/tally;
                err_y = sqrt(total_y2/tally - avg_y*avg_y);
            else
                avg_x = 0;
                avg_y = 0;
                err_y = 0;
            end
            output(counter,:) = [avg_x avg_y err_y];
            counter = counter + 1;
        end
        
        % Output file with a .avg suffix
        out_file = [filename output_extension];
        MiscFns.write_csv(out_file, MiscFns.csv_struct(output, colheaders), delimiter);
        %dlmwrite(out_file, output, '\t');
        disp(sprintf('Processed... Output file: %s', out_file));
        %dlmwrite([filename '.test'], [file.data(:, 2) file.data(:, 3)], '\t');

        % END
    end
end