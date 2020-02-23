clear all
close all
clc


% Rectangular problem area
L = 2.5/100;
H = 2.5/100;

% Grid resolution
resolution = L/100;

delta_x = resolution;
delta_y = resolution;


V_a = 500;
V_b = 100;
V_c = 300;
V_d = 400;



Va = struct(...
    'constant_coordinate', 'y', ...
    'value_of_constant_coordinate', H, ...
    'line_start', 0, ...
    'line_end', L, ...
    'boundary_value', V_a);

Vb = struct(...
    'constant_coordinate', 'x', ...
    'value_of_constant_coordinate', L, ...
    'line_start', delta_x, ...
    'line_end', H-delta_x, ...
    'boundary_value', V_b);

Vc = struct(...
    'constant_coordinate', 'y', ...
    'value_of_constant_coordinate', 0, ...
    'line_start', 0, ...
    'line_end', L, ...
    'boundary_value', V_c);

Vd = struct(...
    'constant_coordinate', 'x', ...
    'value_of_constant_coordinate', 0, ...
    'line_start', delta_x, ...
    'line_end', H-delta_x, ...
    'boundary_value', V_d);

boundary_conditions = {Va,Vb,Vc,Vd};



grid_data = get_FDM_grid(L,H,resolution);
nodes_on_boundary_LUT = specify_boundary(L, H,resolution,boundary_conditions, grid_data);

node_data = FDM(grid_data, nodes_on_boundary_LUT);
%%
postprocessor(L,H,resolution,grid_data, node_data, 1,18);
%%
dir_out = 'Bsp_1_c';
[~,~] = mkdir(dir_out);
h =  findobj('type','figure');
for k = 1 : length(h)
    f = h(k);
    name = fullfile(dir_out,sprintf('fig_%d',f.Number));
    %savefig(h(k),name);
    saveas(h(k),name,'epsc')
    
end
