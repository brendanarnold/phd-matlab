classdef MiscFns
    % A scoped class to hold miscellaneous functions and values
    
    properties (Static = true)
        tab = char(9);
        newline_dos = [char(13) char(10)];
        newline_unix = char(10);
    end
    
    methods (Static = true)        
        
        function [out] = csv_struct(data, headers)
            % A 'constructor' for the data struct used in import_csv and
            % write_csv
            %out = struct('data', [], 'textdata', {}, 'colheaders', {});
            out.data = data;
            out.colheaders = headers;
            out.textdata = {};
        end
        
        function [out] = import_csv(filename, delimiter, header_row)
            % Wrapper around the woeful importdata function which
            % arbitrarily returns a struct or the data depending on the
            % arguments and contents of the file
            if header_row > 0
                % If a header row is called for, importdata returns a
                % struct
                out = importdata(filename, delimiter, header_row);                
            else
                file_data = importdata(filename, delimiter);
                if isstruct(file_data)
                    % Even if header lines are not asked for, importdata return
                    % a struct with headers if a line of strings is found
                    file_data.colheaders = {};
                    file_data.textdata = {};
                    out = file_data;
                else
                    % importdata does not return a struct if the file
                    % contains just numbers and you do not specify headers
                    out = MiscFns.csv_struct(file_data, {});
                end
            end
        end
        
        function write_csv(filename, data, delimiter)
            % Writes files taking the struct defined in 'import_csv' as
            % input. i.e. allows for headers to be included.
            file_handle = fopen(filename, 'w');
            % Write headers to file if there are headers in the struct
            if ~isempty(data.colheaders)
                file_line = '';
                for col_index = 1:length(data.colheaders)
                    if col_index == 1
                        file_line = data.colheaders{col_index};
                    else
                        file_line = [file_line delimiter data.colheaders{col_index}];
                    end
                end
                fprintf(file_handle, '%s', [file_line MiscFns.newline_dos]);
            end
            % Write out the data
            [rows, cols] = size(data.data);
            for row_index = 1:rows
                file_line = '';
                for col_index = 1:cols
                    str_val = num2str(data.data(row_index, col_index));
                    if col_index == 1
                        file_line = str_val;
                    else
                        file_line = [file_line delimiter str_val];
                    end
                end
                fprintf(file_handle, '%s', [file_line MiscFns.newline_dos]);
            end
            % Close the file and actually write the data
            result = fclose(file_handle);
            if result == -1
                error('Problem writing the file')
            end
        end
        
        
        function [out] = ends_with(haystack_string, needle_string)
            % Checks whether the last part of a string matches a substring
            out = 0; % Not matched until proven matched
            if length(haystack_string) >= length(needle_string) % Make sure is long enough to match
                if strcmp(haystack_string(end-length(needle_string)+1:end), needle_string);
                    out = 1;
                end
            end
        end
        
        function [out] = get_cols(matrix, cols)
            % Returns a new matrix containing only the columns specified
            out = [];
            for col_index = cols
                out = [out matrix(:,col_index)];
            end
        end
        
        function [out] = get_files(directory_str, regex_format)
            % Returns all cell array of file strings in given directory that match the regular
            % expression
            out = {};
            file_info = dir(directory_str);
            for index = 1:length(file_info)
                filename = file_info(index).name;
                if MiscFns.regex_match(filename, regex_format)
                    out = [out filename];
                end
            end
        end
        
        function [out] = regex_match(test_str, regex_format)
            % Wrapper around regexp to return a TRUE or FALSE (0 or 1) if
            % the regular expression matches at all
            if isempty(regexp(test_str, regex_format, 'once'))
                out = 0;
            else
                out = 1;
            end
        end
        
        
        
    end
    
end