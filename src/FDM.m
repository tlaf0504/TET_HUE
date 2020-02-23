function node_data = FDM(grid_data, nodes_on_boundary_LUT)

[number_of_nodes,~] = size(grid_data);
[number_of_nodes_on_boundary,~] = size(nodes_on_boundary_LUT);
number_of_unknowns = number_of_nodes - number_of_nodes_on_boundary;


node_number_to_unknown_vector_idx_LUT = zeros(number_of_unknowns, 2);
node_number_to_unknown_vector_idx_LUT(:,2) = 1:number_of_unknowns;

unknown_potential_idx = 1;
for k = 1 : number_of_nodes
    if ~any(nodes_on_boundary_LUT(:,3) == k-1)
        node_number_to_unknown_vector_idx_LUT(unknown_potential_idx,1) = k-1;
        unknown_potential_idx = unknown_potential_idx + 1;
    end
end

FDM_matrix = zeros(number_of_unknowns, number_of_unknowns);
FDM_right_side = zeros(number_of_unknowns, 1);

fprintf('Assembling FDM equation system...')
tic
for k = 1:number_of_unknowns
    idx = node_number_to_unknown_vector_idx_LUT(:,2) == k;
    node_number_of_current_unknown = node_number_to_unknown_vector_idx_LUT(idx,1);
    
    [FDM_matrix_row, FDM_right_side_element] = handle_FDM_element(k, ...
    grid_data, ...
    nodes_on_boundary_LUT, ...
    node_number_to_unknown_vector_idx_LUT);

    FDM_matrix(k,:) = FDM_matrix_row;
    FDM_right_side(k) = FDM_right_side_element;
   
end

fprintf('Done.\nTook %f seconds.\n\n', toc)

fprintf('Solving equation system...')
tic
FDM_solution = FDM_matrix\FDM_right_side;
fprintf('Done.\nTook %f seconds.\n\n', toc)

node_data = get_node_data(grid_data,  ...
    FDM_solution, ...
    nodes_on_boundary_LUT, ...
    node_number_to_unknown_vector_idx_LUT);


end

function [FDM_matrix_row, FDM_right_side] = handle_FDM_element(unknown_number, ...
    grid_data, ...
    nodes_on_boundary, ...
    node_number_to_unknown_vector_idx_LUT)

    [number_of_unknowns,~] =  size(node_number_to_unknown_vector_idx_LUT);

    center_node_number = node_number_to_unknown_vector_idx_LUT(node_number_to_unknown_vector_idx_LUT(:,2)==unknown_number, 1);
    
    tmp = grid_data(grid_data(:,3) == center_node_number, 1:2);
    row_idx = tmp(1);
    col_idx = tmp(2);
    
    max_row_idx = max(grid_data(:,1));
    max_col_idx = max(grid_data(:,2));

    FDM_matrix_row = zeros(1, number_of_unknowns);
    FDM_right_side = 0;
    
    % Center element
    FDM_matrix_row(unknown_number) = 4;
    
    % Upper element
    current_row_idx = row_idx-1;
    current_col_idx = col_idx;
    boundary_value = check_if_node_is_on_boundary(current_row_idx, current_col_idx, nodes_on_boundary);
    
    if ~isnan(boundary_value)
        FDM_right_side = FDM_right_side + boundary_value;
    else
        if current_row_idx >= 0 && current_col_idx >= 0 && ...
           current_row_idx <= max_row_idx && current_col_idx <= max_col_idx     
            idx = (grid_data(:,1) == current_row_idx) & (grid_data(:,2) == current_col_idx);
            current_node_number = grid_data(idx,3);
            current_node_unknown_number = node_number_to_unknown_vector_idx_LUT(node_number_to_unknown_vector_idx_LUT(:,1) == current_node_number, 2);
            FDM_matrix_row(current_node_unknown_number) = -1;
        end
    end
    
    % Lower element
    current_row_idx = row_idx+1;
    current_col_idx = col_idx;
    boundary_value = check_if_node_is_on_boundary(current_row_idx, current_col_idx, nodes_on_boundary);
    
    if ~isnan(boundary_value)
        FDM_right_side = FDM_right_side + boundary_value;
    else
        if current_row_idx >= 0 && current_col_idx >= 0 && ...
           current_row_idx <= max_row_idx && current_col_idx <= max_col_idx     
            idx = (grid_data(:,1) == current_row_idx) & (grid_data(:,2) == current_col_idx);
            current_node_number = grid_data(idx,3);
            current_node_unknown_number = node_number_to_unknown_vector_idx_LUT(node_number_to_unknown_vector_idx_LUT(:,1) == current_node_number, 2);
            FDM_matrix_row(current_node_unknown_number) = -1;
        end
    end
    
    
    
    % Left element
    current_row_idx = row_idx;
    current_col_idx = col_idx-1;
    boundary_value = check_if_node_is_on_boundary(current_row_idx, current_col_idx, nodes_on_boundary);
    
    if ~isnan(boundary_value)
        FDM_right_side = FDM_right_side + boundary_value;
    else
        if current_row_idx >= 0 && current_col_idx >= 0 && ...
           current_row_idx <= max_row_idx && current_col_idx <= max_col_idx     
            idx = (grid_data(:,1) == current_row_idx) & (grid_data(:,2) == current_col_idx);
            current_node_number = grid_data(idx,3);
            current_node_unknown_number = node_number_to_unknown_vector_idx_LUT(node_number_to_unknown_vector_idx_LUT(:,1) == current_node_number, 2);
            FDM_matrix_row(current_node_unknown_number) = -1;
        end
    end
    
    
    % Right element
    current_row_idx = row_idx;
    current_col_idx = col_idx+1;
    boundary_value = check_if_node_is_on_boundary(current_row_idx, current_col_idx, nodes_on_boundary);
    
    if ~isnan(boundary_value)
        FDM_right_side = FDM_right_side + boundary_value;
    else
        if current_row_idx >= 0 && current_col_idx >= 0 && ...
           current_row_idx <= max_row_idx && current_col_idx <= max_col_idx     
            idx = (grid_data(:,1) == current_row_idx) & (grid_data(:,2) == current_col_idx);
            current_node_number = grid_data(idx,3);
            current_node_unknown_number = node_number_to_unknown_vector_idx_LUT(node_number_to_unknown_vector_idx_LUT(:,1) == current_node_number, 2);
            FDM_matrix_row(current_node_unknown_number) = -1;
        end
    end

end

function boundary_value = check_if_node_is_on_boundary(row_idx, col_idx, nodes_on_boundary)

    if any((nodes_on_boundary(:,1) == row_idx) & (nodes_on_boundary(:,2) == col_idx))
        idx = (nodes_on_boundary(:,1) == row_idx) & (nodes_on_boundary(:,2) == col_idx);
        boundary_value = nodes_on_boundary(idx, 4);
        
    else
        boundary_value = nan;
    end

end
    
function node_data = get_node_data(grid_data,  ...
    FDM_solution, ...
    nodes_on_boundary_LUT, ...
    node_number_to_unknown_vector_idx_LUT)

    [number_of_nodes,~] = size(grid_data);
    
    node_data = zeros(number_of_nodes, 1);
    
    for k = 0 : number_of_nodes-1
        if any(nodes_on_boundary_LUT(:,3) == k)
            idx = nodes_on_boundary_LUT(:,3) == k;
            node_data(k+1) = nodes_on_boundary_LUT(idx,4);
        else
            idx = node_number_to_unknown_vector_idx_LUT(:, 1) == k;
            idx_2 = node_number_to_unknown_vector_idx_LUT(idx, 2);
            node_data(k+1) = FDM_solution(idx_2);
        end
    end
end


