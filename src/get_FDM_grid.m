function row_col_idx_to_node_number_LUT = get_FDM_grid(L,H,resolution)

steps_x = 0:resolution:L;
steps_y = 0:resolution:H;

N_steps_x = length(steps_x);
N_steps_y = length(steps_y);

global_node_numbers = 0:(N_steps_x*N_steps_y)-1;

colum_numbers = 0:(N_steps_x)-1;

row_col_idx_to_node_number_LUT = zeros((N_steps_x*N_steps_y), 3);

for k = 1 : N_steps_y
    idx_range = (k-1)*(N_steps_x) + (1:N_steps_x);
    row_col_idx_to_node_number_LUT(idx_range, 1) = (k-1);
    row_col_idx_to_node_number_LUT(idx_range,2) = colum_numbers;
end

row_col_idx_to_node_number_LUT(:,3) = global_node_numbers;

end