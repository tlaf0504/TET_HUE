function nodes_on_boundary_LUT = specify_boundary(L,H,resolution,boundary_conditions, grid_data)

number_of_boundary_conditions = length(boundary_conditions);

delta_x = resolution;
delta_y = resolution;

nodes_on_boundary_LUT = [];

for k = 1 : number_of_boundary_conditions
    current_boundary_condition = boundary_conditions{k};
    
    % Line in y-direction
    if strcmpi(current_boundary_condition.constant_coordinate, 'x')
        col_idx_of_line = round(...
            current_boundary_condition.value_of_constant_coordinate/delta_x);
        
        start_row_idx_of_line = round(...
            current_boundary_condition.line_start/delta_y);
        end_row_idx_of_line = round(...
            current_boundary_condition.line_end/delta_y);
        
        
        row_idx_range = start_row_idx_of_line : end_row_idx_of_line;
        
        for l = 1:length(row_idx_range)
            logical_row_idx = (grid_data(:,1) ==  row_idx_range(l)) & ...
                (grid_data(:,2) == col_idx_of_line);
            
            node_number = grid_data(logical_row_idx,3);
            nodes_on_boundary_LUT = [nodes_on_boundary_LUT;...
                grid_data(logical_row_idx,1), grid_data(logical_row_idx,2),node_number, current_boundary_condition.boundary_value];
        end
        
    % Line in x-direction    
    elseif strcmpi(current_boundary_condition.constant_coordinate, 'y')
        
        row_idx_of_line = round(...
            current_boundary_condition.value_of_constant_coordinate/delta_y);
        
        start_col_idx_of_line = round(...
            current_boundary_condition.line_start/delta_x);
        end_col_idx_of_line = round(...
            current_boundary_condition.line_end/delta_x);
        
        
        col_idx_range = start_col_idx_of_line : end_col_idx_of_line;
        
        for l = 1:length(col_idx_range)
            logical_row_idx = (grid_data(:,1) == row_idx_of_line) & ...
                (grid_data(:,2) == col_idx_range(l));
            
            node_number = grid_data(logical_row_idx,3);
            nodes_on_boundary_LUT = [nodes_on_boundary_LUT;...
                grid_data(logical_row_idx,1), grid_data(logical_row_idx,2),node_number, current_boundary_condition.boundary_value];
        end
        
        
    else
        error(['Invalid boundary condition specification. Use either "x"', ...
            'or "y" for specifying constant coordinate'])
    end
end


% Check for nodes on multiple boundaries
node_numbers = sortrows(nodes_on_boundary_LUT(:,3));
d = diff(node_numbers);
if any(d == 0)
    error('Some nodes are defined on multiple boundaries. Check definition of boundary conditions.')
end

end

