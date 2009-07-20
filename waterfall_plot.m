config_filename = 'waterfall_config.txt';
input_extension = '.*_(.*)-.*_spliced.dat$';

if ~exist(config_filename, 'file')
    dir_output = dir();
    file_info = {};
    file_line = 1;
    for index = 1:length(dir_output)
        filename = dir_output(index).name; % Cache the filename
        [legend] = regexp(filename, input_extension, 'tokens');
        if ~isempty(legend) % If the filename matched the pattern ...
            file_info{file_line, 1} = filename;
            file_info{file_line, 2} = legend{1};
            file_line = file_line + 1;
        end
    end
    file_data = MiscFns.csv_struct(file_info, {'filename' 'legend'});
    MiscFns.write_csv(config_filename, file_data, MiscFns.tab);
end
edit(config_filename);
display('Push any key once finished editing the file configuration ...');
pause
% Wait around for the user to edit the config file...

display('Plot the thing!');

file_handle = fopen(config_filename);
cell_array = textscan(file_handle, '%s %s', 'Delimiter', MiscFns.tab);


