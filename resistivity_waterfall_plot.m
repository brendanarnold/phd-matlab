function resistivity_waterfall_plot(y_col)
    % Generate a 2D waterfall plot from an automatically generated
    % configuration file.
    config_filename = 'waterfall_config.txt';
    input_extension = '.*_(.*)-.*_spliced\.dat\.avg$';
%    input_extension = '.*_(.*)\.hall$';
%    config_filename = 'hall_waterfall_config.txt';
    x_col = 1;
    %y_col = 4;
    
    if ~exist(config_filename, 'file')
        dir_output = dir();
        file_info = {};
        file_line = 1; % First line is for the title
        for index = 1:length(dir_output)
            filename = dir_output(index).name; % Cache the filename
            [legend_str] = regexp(filename, input_extension, 'tokens');
            if ~isempty(legend_str) % If the filename matched the pattern ...
                file_info{file_line, 1} = filename;
                file_info{file_line, 2} = legend_str{1};
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
    line_colors = [];
    num_temps = 19;
    for index = 1:length(legends)
        switch lower(legends{index})
            case '1p5k'
                line_colors(index,:) = MiscFns.extract_color(1/num_temps);
            case '4p2k'
                line_colors(index,:) = MiscFns.extract_color(2/num_temps);
            case '10k'
                line_colors(index,:) = MiscFns.extract_color(3/num_temps);
            case '15k'
                line_colors(index,:) = MiscFns.extract_color(4/num_temps);
            case '20k'
                line_colors(index,:) = MiscFns.extract_color(4/num_temps);
            case '25k'
                line_colors(index,:) = MiscFns.extract_color(5/num_temps);
            case '30k'
                line_colors(index,:) = MiscFns.extract_color(6/num_temps);
            case '35k'
                line_colors(index,:) = MiscFns.extract_color(7/num_temps);
            case '40k'
                line_colors(index,:) = MiscFns.extract_color(8/num_temps);
            case '45k'
                line_colors(index,:) = MiscFns.extract_color(9/num_temps);
            case '50k'
                line_colors(index,:) = MiscFns.extract_color(10/num_temps);
            case '55k'
                line_colors(index,:) = MiscFns.extract_color(11/num_temps);
            case '60k'
                line_colors(index,:) = MiscFns.extract_color(12/num_temps);
            case '70k'
                line_colors(index,:) = MiscFns.extract_color(13/num_temps);
            case '75k'
                line_colors(index,:) = MiscFns.extract_color(14/num_temps);
            case '80k'
                line_colors(index,:) = MiscFns.extract_color(15/num_temps);
            case '100k'
                line_colors(index,:) = MiscFns.extract_color(16/num_temps);
            case '120k'
                line_colors(index,:) = MiscFns.extract_color(17/num_temps);
            case '140k'
                line_colors(index,:) = MiscFns.extract_color(18/num_temps);
            case '160k'
                line_colors(index,:) = MiscFns.extract_color(19/num_temps);
        end
    end


    % Plot line
    num_plots = length(filenames);
    lines = zeros(1, num_plots);
    figure();
    for index = 1:num_plots
        %display(filenames{index});
        file_struct = MiscFns.import_csv(filenames{index}, MiscFns.tab, 1);
        line_color = line_colors(index,:);
        lines(index) = line(file_struct.data(:,x_col), file_struct.data(:,y_col), 'Color', line_color);
        display(sprintf('Drawing %s ...', filenames{index}));
    end
    legend(lines, legends, 'Position', [0.6561 0.1929 0.1533 0.3333]);
    xlabel('B (T)');
    ylabel('Rhall');
    title(plot_title);

    % Some feedback
    display(sprintf('2D waterfall plot `%s` generated!', plot_title));

end