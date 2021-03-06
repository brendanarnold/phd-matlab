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
                    val = data.data(row_index, col_index);
                    if isnumeric(val) % Convert numeric to string
                        str_val = num2str(val);
                    elseif iscell(val) % Convert cell to string
                        while iscell(val) % Why do I need to do this?
                            val = val{1};
                        end
                        str_val = val; 
                    else % Assume all else is string (not true I know)
                        str_val = val;
                    end
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
        
%         function [out] = regex_match(test_str, regex_format)
%             % Wrapper around regexp to return a TRUE or FALSE (0 or 1) if
%             % the regular expression matches at all
%             if isempty(regexp(test_str, regex_format, 'once'))
%                 out = 0;
%             else
%                 out = 1;
%             end
%         end
%         
%         %function serialize_structs(struct_array, out_filename)
%         %    MiscFns.serialize_structs(struct_array, 'root', 'datafile', out_filename);
%         %end
%         
%         function serialize_structs(struct_array, root_element_name, struct_element_name, out_filename)
%             % Serializes an array of structs into an xml document
%             % This is useful for creating a configuration file
%             % Caveats:
%             %  Can only deal with fields containing strings
%             %  Can only serialise structs one level deep
%             document_handle = com.mathworks.xml.XMLUtils.createDocument(root_element_name);
%             root_node = document_handle.getDocumentElement;
%             for index = 1:length(struct_array)
%                 struct_element = document_handle.createElement(struct_element_name);
%                 for field_name = fieldnames(struct_array(index))'
%                     field_name = field_name{1};
%                     field_str = getfield(struct_array(index), field_name);
%                     field_element = document_handle.createElement(field_name);
%                     field_element.appendChild(document_handle.createTextNode(field_str));
%                     struct_element.appendChild(field_element);
%                 end
%                 root_node.appendChild(struct_element);
%             end
%             xmlwrite(out_filename, document_handle);
%         end
%         
%         function [out] = unserialize_structs(filename)
%             document_handle = xmlread(filename);
%             document_handle.
%         end
        
        function [out] = extract_color(frctn, cmap)
            % Extract closest color from default colormap or given colormap given
            % a floating point number between 0 and 1.
            if ~exist('cmap', 'var') || ~isnumeric(cmap) || ~isempty(cmap) 
                cmap = colormap();
            end
            index = ceil(frctn * size(cmap, 1)); % Ceil because matlab is not zero indexed
            if index == 0
                out = cmap(1,:);
            else
                out = cmap(index,:);
            end
        end
        
        
        
    end
    
end