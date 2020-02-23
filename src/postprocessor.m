function postprocessor(L,H,resolution,grid_data, node_data,material,N_skip)


boundary_lines_x = [0,L,L,0,0];
boundary_lines_y = [0,0,H,H,0];


fprintf('Starting postprocessor...')
tic
delta_x = resolution;
delta_y = resolution;

N_steps_x = L/resolution + 1;
N_steps_y = H/resolution + 1;

x = 0:delta_x:L;
y = 0:delta_y:H;

[X,Y] = meshgrid(x,y);
Z = reshape(node_data, N_steps_x, N_steps_y)';

figure('Name', 'Potential')
surf(X,Y,Z,'EdgeColor', 'none');
grid minor
xlabel('Distance in m')
ylabel('Distance in m')
zlabel('Potential in V')
title('Potential')
colormap jet
colorbar

figure('Name', 'Potential with equipotential lines')
surf(X,Y,Z, 'EdgeColor', 'none', 'FaceAlpha', 1);
hold on
contour(X,Y,Z,'LineWidth',0.5, 'Linecolor', 'black')
view(180,-90);
xlabel('Distance in m')
ylabel('Distance in m')
zlabel('Potential in V')
title('Potential with equipotential lines')
colormap jet
colorbar


% =========================== Calculate field data ========================

field_data = zeros(N_steps_x-1*N_steps_y-1, 1);
field_data_2 = zeros(N_steps_x-1, N_steps_y-1);
field_coordinates = zeros(N_steps_x-1*N_steps_y-1, 1);

for k = 0:N_steps_y-2
    for l = 0 : N_steps_x-2
        idx = (k*N_steps_x) + l + 1;
        field_coordinates(idx,1) = l*resolution + resolution/2;
        field_coordinates(idx,2) = k*resolution + resolution/2;
        
        
        tmp = get_field_at_point_xy(l*resolution, k*resolution,grid_data, node_data, material, resolution);
        field_data(idx, 1) = tmp(1);
        field_data(idx, 2) = tmp(2);
        
        field_data_2(l+1,k+1) = tmp(1) + 1i*tmp(2);
        
    end
end


% =========================== Plot field ==================================
abs_field = sqrt(field_data(:,1).^2 + field_data(:,2).^2);
scale = 1/max(abs_field)*L/10;
field_data = field_data*scale;

field_dir_u = field_data(:,1)./abs_field;
field_dir_v = field_data(:,2)./abs_field;

field_x = field_coordinates(1:N_skip:end,1); % Don't use all field arrows
field_y = field_coordinates(1:N_skip:end,2);
field_u = field_dir_u(1:N_skip:end);
field_v = field_dir_v(1:N_skip:end);

figure('Name', 'Field direction')
quiver(field_x, field_y, field_u,field_v)
hold on
plot(boundary_lines_x, boundary_lines_y, 'Color', 'black', 'LineWidth', 2)
axis equal
grid minor
title('Field direction')
xlabel('x in m')
ylabel('y in m')

% figure('Name', 'Field direction')
% quiverwcolorbar(field_x, field_y, field_u, field_v)
% hold on
% plot(boundary_lines_x, boundary_lines_y, 'Color', 'black', 'LineWidth', 2)
% axis equal
% grid minor
% title('Field direction')
% xlabel('x in m')
% ylabel('y in m')
% 
% 
% 
% 
% figure('Name', 'Field direction')
% quiverC2D(field_x, field_y, field_u, field_v)
% hold on
% plot(boundary_lines_x, boundary_lines_y, 'Color', 'black', 'LineWidth', 2)
% axis equal
% grid minor
% title('Field direction')
% xlabel('x in m')
% ylabel('y in m')




X = X'; Y = Y';
Z_field = griddata(X(1:end-1,1:end-1),Y(1:end-1,1:end-1),abs(field_data_2),X,Y);
%Z_field = abs(field_data_2);

figure('Name', 'Field intensity')
surf(X,Y,Z_field,'EdgeColor','none');
grid minor
title('Field intensity')
xlabel('x in m')
ylabel('y in m')
zlabel('Field strength')
colormap jet
colorbar

fprintf('Done.\nTook %f seconds.\n', toc)




end

function V = get_potential_at_point_xy(x,y,grid_data,node_data, resolution)

[boundary_node_coordinates, boundary_node_potentials] = ...
    get_current_element_information(x,y,grid_data, node_data,resolution);

V1 = boundary_node_potentials(1);
V2 = boundary_node_potentials(2);
V3 = boundary_node_potentials(3);
V4 = boundary_node_potentials(4);


tmp = ([x,y] - boundary_node_coordinates(1,:))./resolution;
xi = tmp(1);
eta = tmp(2);


V = V1*(1-xi)*(1-eta) + V2*(1-xi)*eta + V3*xi*eta + V4*(1-xi)*eta;

end

function E = get_field_at_point_xy(x,y,grid_data, node_data, material,resolution)
E = -material*grad_V(x,y,grid_data,node_data,resolution);


end

function vec = grad_V(x,y,grid_data,node_data, resolution)
[boundary_node_coordinates, boundary_node_potentials] = ...
    get_current_element_information(x,y,grid_data, node_data,resolution);

V1 = boundary_node_potentials(1);
V2 = boundary_node_potentials(2);
V3 = boundary_node_potentials(3);
V4 = boundary_node_potentials(4);


tmp = ([x,y] - boundary_node_coordinates(1,:))./resolution;
xi = tmp(1);
eta = tmp(2);

dN1_dXi = -(1-eta);
dN1_dEta = -(1-xi);

dN2_dXi = -eta;
dN2_dEta = (1-xi);

dN3_dXi = eta;
dN3_dEta = xi;

dN4_dXi = (1-eta);
dN4_dEta = -xi;

dXi_dX = 1/resolution;
dXi_dY = 0;

dEta_dX = 0;
dEta_dY = 1/resolution;

dV_dXi = V1*dN1_dXi + V2*dN2_dXi + V3*dN3_dXi + V4*dN4_dXi;
dV_dEta = V1*dN1_dEta + V2*dN2_dEta + V3*dN3_dEta + V4*dN4_dEta;

vec = [ ...
    dV_dXi*dXi_dX + dV_dEta*dEta_dX; ...
    dV_dXi*dXi_dY + dV_dEta*dEta_dY];

end


function [boundary_node_coordinates, boundary_node_potentials] = ...
    get_current_element_information(x,y,grid_data,node_data,resolution)


upper_left_node_relative_x_coordinate = floor(x/resolution);
upper_left_node_relative_y_coordinate = floor(y/resolution);

lower_left_node_relative_x_coordinate = upper_left_node_relative_x_coordinate;
lower_left_node_relative_y_coordinate = upper_left_node_relative_y_coordinate+1;

upper_right_node_relative_x_coordinate = upper_left_node_relative_x_coordinate+1;
upper_right_node_relative_y_coordinate = upper_left_node_relative_y_coordinate;

lower_right_node_relative_x_coordinate = upper_left_node_relative_x_coordinate+1;
lower_right_node_relative_y_coordinate = upper_left_node_relative_y_coordinate+1;

upper_left_node_potential_idx = ...
    (grid_data(:,2) == upper_left_node_relative_x_coordinate) & ...
    (grid_data(:,1) == upper_left_node_relative_y_coordinate);

lower_left_node_potential_idx = ...
    (grid_data(:,2) == lower_left_node_relative_x_coordinate) & ...
    (grid_data(:,1) == lower_left_node_relative_y_coordinate);

upper_right_node_potential_idx = ...
    (grid_data(:,2) == upper_right_node_relative_x_coordinate) & ...
    (grid_data(:,1) == upper_right_node_relative_y_coordinate);

lower_right_node_potential_idx = ...
    (grid_data(:,2) == lower_right_node_relative_x_coordinate) & ...
    (grid_data(:,1) == lower_right_node_relative_y_coordinate);



upper_left_node_potential = node_data(upper_left_node_potential_idx);
lower_left_node_potential = node_data(lower_left_node_potential_idx);
upper_right_node_potential = node_data(upper_right_node_potential_idx);
lower_right_node_potential = node_data(lower_right_node_potential_idx);

upper_left_absolute_node_coordinates = ...
    [upper_left_node_relative_x_coordinate,upper_left_node_relative_y_coordinate] * ...
    resolution;

lower_left_absolute_node_coordinates = ...
    [lower_left_node_relative_x_coordinate,lower_left_node_relative_y_coordinate] * ...
    resolution;

upper_right_absolute_node_coordinates = ...
    [upper_right_node_relative_x_coordinate,upper_right_node_relative_y_coordinate] * ...
    resolution;

lower_right_absolute_node_coordinates = ...
    [lower_right_node_relative_x_coordinate,lower_right_node_relative_y_coordinate] * ...
    resolution;

boundary_node_coordinates = [ ...
    upper_left_absolute_node_coordinates; ...
    lower_left_absolute_node_coordinates; ...
    lower_right_absolute_node_coordinates; ...
    upper_right_absolute_node_coordinates];


boundary_node_potentials = [ ...
    upper_left_node_potential; ...
    lower_left_node_potential; ...
    lower_right_node_potential; ...
    upper_right_node_potential
    ];

end

function positions = get_equipotential_lines(grid_data, node_data, resolution, L, H)

N_points=10;
y_start = linspace(resolution,H-resolution,N_points);
x_start = linspace(resolution,L-resolution,N_points);
positions = cell(N_points, 1);

for k = 0 : N_points-1
     positions{k+1} = get_single_equipotential_line(0,y_start(k+1), [1;0], L, H, grid_data, node_data, resolution);
     
end

% for k = 0 : N_points-1
%      positions{N_points + k+1} = ...
%          get_single_equipotential_line(x_start(k+1),H-resolution, [0;-1], L, H, grid_data, node_data, resolution);
%      
% end

% for k = 0 : N_points-1
%      positions{2*N_points + k+1} = ...
%          get_single_equipotential_line(L-resolution, y_start(k+1), [-1;0], L, H, grid_data, node_data, resolution);
%      
% end,potential_values,
% 
% for k = 0 : N_points-1
%      positions{3*N_points + k+1} = ...
%          get_single_equipotential_line(x_start(k+1),resolution, [0;-1], L, H, grid_data, node_data, resolution);
%      
% end


end

function positions = get_single_equipotential_line(x_start, y_start, e_start, L, H, grid_data, node_data, resolution)

x_cur = x_start;
y_cur = y_start;

positions = [x_cur, y_cur];

eim1 = e_start;


while (x_cur <= L) && (y_cur <= H) && (x_cur >= 0) && (y_cur >= 0)
    
    vec = grad_V(x_cur,y_cur,grid_data,node_data, resolution);
    
    Ei = -vec/norm(vec);
    Eix = Ei(1);
    Eiy = Ei(2);
    
    if Eix ~= 0 && Eiy ~= 0
        eix = -Eiy/Eix;
        ei = [eix; 1];
        ei = ei/norm(ei);
    elseif Eix == 0 && Eiy~= 0
        ei = [1;0];
    elseif Eix ~= 0 && Eiy == 0
        ei = [0;1];
    else
        ei = [0;0];
        error('Field at point is [0,0]. Could nor determine euquipotential line')
    end
    
    if ei'*eim1 < 0
        ei = -ei;
    end
    
    positions = [positions; [x_cur, y_cur]];
    
    x_cur = x_cur + ei(1)*resolution;
    y_cur = y_cur + ei(2)*resolution;
    
%     fprintf('x=%f\ty=%f\n', x_cur, y_cur)
%     disp(Ei)
%     disp(ei)
%     fprintf('\n\n')
    eim1 = ei;
    
    
end

end


