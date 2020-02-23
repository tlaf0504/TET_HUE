clear all
close all
clc


% Rectangular problem area
L = 12e-3;
H = 5e-3;

% Grid resolution
resolution = 1e-4;


delta_x = resolution;
delta_y = resolution;

Va = struct(...
    'constant_coordinate', 'y', ...
    'value_of_constant_coordinate', 0 , ...
    'line_start', 0, ...
    'line_end', L, ...
    'boundary_value', 0);



Vb = struct(...
    'constant_coordinate', 'y', ...
    'value_of_constant_coordinate', H, ...
    'line_start', 4e-3, ...
    'line_end', 8e-3, ...
    'boundary_value', 10);




boundary_conditions = {Va,Vb};


grid_data = get_FDM_grid(L,H,resolution);
nodes_on_boundary_LUT = specify_boundary(L, H,resolution,boundary_conditions, grid_data);

node_data = FDM(grid_data, nodes_on_boundary_LUT);
postprocessor(L,H,resolution,grid_data, node_data, 3.7e7,16);

%%
dir_out = 'Bsp_2';
[~,~] = mkdir(dir_out);
h =  findobj('type','figure');
for k = 1 : length(h)
    f = h(k);
    name = fullfile(dir_out,sprintf('fig_%d',f.Number));
    %savefig(h(k),name);
    saveas(h(k),name,'epsc')
    
end