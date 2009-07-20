config_filename = 'waterfall_config.txt';
input_extension = '.*_(.*)-.*_spliced.dat.avg$';

if ~exist(config_filename, 'file')
    dir_output = dir();
    file_info = {};
    file_line = 1; % First line is for the title
    for index = 1:length(dir_output)
        filename = dir_output(index).name; % Cache the filename
        [legend] = regexp(filename, input_extension, 'tokens');
        if ~isempty(legend) % If the filename matched the pattern ...
            file_info{file_line, 1} = filename;
            file_info{file_line, 2} = legend{1};
            file_line = file_line + 1;
        end
    end
    file_data = MiscFns.csv_struct(file_info, {'title:' 'TITLE HERE'}); % Not really the title but used for convenience
    MiscFns.write_csv(config_filename, file_data, MiscFns.tab);
end
edit(config_filename);
display('Push any key once finished editing the file configuration ...');
pause
% Wait around for the user to edit the config file...

display('Plot the thing ...');

% Read in what the configuration file says
file_handle = fopen(config_filename);
file_info = textscan(file_handle, '%s %s', 'Delimiter', MiscFns.tab); 
% This returns a 1x2 cell array, each cell contaning another nx1 cell
% array. Why?!
fclose(file_handle);
first_column = file_info{1};
filenames = first_column(2:end);
second_column = file_info{2};
plot_title = second_column{1}; 
legends = second_column(2:end);


% Plot line
num_plots = length(filenames);
legend_cell = {};
lines = zeros(1, num_plots);
figure();
for index = 1:num_plots
    %display(filenames{index});
    file_struct = MiscFns.import_csv(filenames{index}, MiscFns.tab, 1);
    color = MiscFns.extract_color(index/num_plots);
    lines(index) = line(file_struct.data(:,1), file_struct.data(:,2), 'Color', color);
    display(sprintf('Drawing %s ...', filenames{index}));
end
legend(lines, legends, 'Position', [0.6561 0.1929 0.1533 0.3333]);
xlabel('B (T)');
ylabel('Resistance (Ohm)');
title(plot_title);

% Some feedback
display(sprintf('2D waterfall plot `%s` generated!', plot_title));
