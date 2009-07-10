classdef ContourFns
    % Use the ContourFns classname to scope functions

    methods (Static = true)

        function contour_output = contourc_plaid(X, Y, Z, height)
            % Wrapper around contourc so as to allow plaid format used
            % for 'contour' to also be used for 'contourc'
            contour_output = contourc(X(1,:), Y(:,1)', Z, [height height]);
        end
        
        function area = calc_area(contour_output)
            % Calculates the output of a captured contour() or
            % contourc() call into an area

            %contour_output = contourc(X, Y, Z, [height height]);
            area = 0.0;
            shapes = ContourFns.convert_to_shapes(contour_output);

            % Cache the length when used in for loop
            len_shapes = length(shapes);
            % Iterating over a cellarray gives back 'cells' so use indexing
            for index = 1:len_shapes
                area = area + polyarea(shapes{index}(1,:), shapes{index}(2,:));
            end
        end
        
        function parsed = convert_to_shapes(contour_output)
            % Converts the output of a captured contour() or contourc() call into a
            % Cell array of matrices. For each matrix, the first row holds x vals,
            % second holds y vals for all vertices in that shape. Each cell
            % contains one shape.

            % Cache the length of the array (referred to in 'for' loop)
            len_contour_output = size(contour_output);
            len_contour_output = len_contour_output(2);
            % Set the current index
            next_index_metadata = 1;
            shape_num = 1;
            index_offset = 1;
            vertices_buffer = [];
            parsed = cell(0);

            for index = 1:len_contour_output
                % Read out the number of vertices in shape from the metadata if the index point is
                % reached in the data
                if index == next_index_metadata
                    next_index_metadata = next_index_metadata + contour_output(2, index) + 1;
                    % Some actions that are not performed when reading the first
                    % set of meta data values
                    if index > 1
                        index_offset = index;
                        parsed{shape_num} = vertices_buffer;
                        shape_num = shape_num + 1;
                        vertices_buffer = [];
                    end
                    % Otherwise store the vertex in the correct array
                else
                    vertices_buffer(:,index-index_offset) = contour_output(:,index);
                end
            end
            % Double make sure right thing is returned
            if ~isempty(vertices_buffer)
                parsed{shape_num} = vertices_buffer;
            end
        end
        
    end

end