clear all
close all
clc

task_var = 'C';
VA = 500;
VB = 100;
VC = 300;
VD = 400;

 L = 2.5/100;
 
 
 N_sum_steps = 110;
 
 N_points = 200;
 tmp = linspace(0,L,N_points);
 
 [X,Y] = meshgrid(tmp,tmp);
 
 V1 = zeros(N_points,N_points);
 V2 = zeros(N_points,N_points);
 V3 = zeros(N_points,N_points);
 V4 = zeros(N_points,N_points);
 
 for k = 1 : N_points
     for l = 1:N_points
         [V1(k,l), V2(k,l), V3(k,l), V4(k,l)] = calc_V (X(k,l),Y(k,l),N_sum_steps,L,VA,VB,VC,VD);
     end
 end

 %% Rotated view
 figure()
 p1 = subplot(221);
 surf(X*100,Y*100,V1,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V1,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_A')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_A in V')
 set(gca,'FontSize',18), set(gcf, 'Position', get(0, 'Screensize'));
 
 p2 = subplot(222);
 surf(X*100,Y*100,V2,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V2,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_B')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_B in V')
 set(gca,'FontSize',18)
 
 p3 = subplot(223);
 surf(X*100,Y*100,V3,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V3,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_C')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_C in V')
 set(gca,'FontSize',18)
 
 p4 = subplot(224);
 surf(X*100,Y*100,V4,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V4,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_D')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_D in V')
 set(gca,'FontSize',18)
 
 figure()
 surf(X*100,Y*100,V1+V2+V3+V4,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V1+V2+V3+V4,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_{ges}')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_{ges} in V')
 set(gca,'FontSize',18), set(gcf, 'Position', get(0, 'Screensize'));
 
 
 %% Standard view
 figure()
 set(gcf, 'Position', get(0, 'Screensize'));
 p1 = subplot(221);
 surf(X*100,Y*100,V1,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V1,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_A')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_A in V')
 set(gca,'FontSize',18)
 view(2)
 
 p2 = subplot(222);
 surf(X*100,Y*100,V2,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V2,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_B')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_B in V')
 set(gca,'FontSize',18)
 view(2)
 
 p3 = subplot(223);
 surf(X*100,Y*100,V3,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V3,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_C')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_C in V')
 set(gca,'FontSize',18)
 view(2)
 
 p4 = subplot(224);
 surf(X*100,Y*100,V4,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V4,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_D')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_D in V')
 set(gca,'FontSize',18)
 view(2)
 
 figure()
 set(gcf, 'Position', get(0, 'Screensize'));
 surf(X*100,Y*100,V1+V2+V3+V4,'EdgeColor','none');
 hold on, contour3(X*100,Y*100,V1+V2+V3+V4,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 colorbar
 colormap jet
 title('V_{ges}')
 xlabel('x in cm'), ylabel('y in cm'), zlabel('V_{ges} in V')
 set(gca,'FontSize',18)
 view(2)
 
 
 
 %% Field plots
 hx = X(1,2)-X(1,1);
 hy = Y(2,1)-Y(1,1);
 
 figure()
 set(gcf, 'Position', get(0, 'Screensize'));
 subplot(221)     
     [gx,gy] = gradient(-V1,hx,hy);
     quiver(100*X(1:size(X,1)/25:end,1:size(X,2)/25:end),...
         100*Y(1:size(Y,1)/25:end,1:size(Y,2)/25:end),...
         gx(1:size(gx,1)/25:end,1:size(gx,2)/25:end),...
         gy(1:size(gy,1)/25:end,1:size(gy,2)/25:end),100);
     hold on, contour(X*100,Y*100,V1,'LineWidth',0.5, 'Linecolor', 'black')
     hold off
     xlim(100*[0,L]), ylim(100*[0,L])
     title('V_A aktiv')
     xlabel('x in cm'), ylabel('y in cm')
     set(gca,'FontSize',18)
     
 subplot(2,2,2)
     [gx,gy] = gradient(-V2,hx,hy);
     quiver(100*X(1:size(X,1)/25:end,1:size(X,2)/25:end),...
         100*Y(1:size(Y,1)/25:end,1:size(Y,2)/25:end),...
         gx(1:size(gx,1)/25:end,1:size(gx,2)/25:end),...
         gy(1:size(gy,1)/25:end,1:size(gy,2)/25:end),25);
     hold on, contour(X*100,Y*100,V2,'LineWidth',0.5, 'Linecolor', 'black')
     hold off
     xlim(100*[0,L]), ylim(100*[0,L])
     title('V_B aktiv')
     xlabel('x in cm'), ylabel('y in cm')
     set(gca,'FontSize',18)
     
  subplot(2,2,3)
     [gx,gy] = gradient(-V3,hx,hy);
     quiver(100*X(1:size(X,1)/25:end,1:size(X,2)/25:end),...
         100*Y(1:size(Y,1)/25:end,1:size(Y,2)/25:end),...
         gx(1:size(gx,1)/25:end,1:size(gx,2)/25:end),...
         gy(1:size(gy,1)/25:end,1:size(gy,2)/25:end),25);
     hold on, contour(X*100,Y*100,V3,'LineWidth',0.5, 'Linecolor', 'black')
     hold off
     xlim(100*[0,L]), ylim(100*[0,L])
     title('V_C aktiv')
     xlabel('x in cm'), ylabel('y in cm')
     set(gca,'FontSize',18)
     
  subplot(2,2,4)
     [gx,gy] = gradient(-V4,hx,hy);
     quiver(100*X(1:size(X,1)/25:end,1:size(X,2)/25:end),...
         100*Y(1:size(Y,1)/25:end,1:size(Y,2)/25:end),...
         gx(1:size(gx,1)/25:end,1:size(gx,2)/25:end),...
         gy(1:size(gy,1)/25:end,1:size(gy,2)/25:end),50);
     hold on, contour(X*100,Y*100,V4,'LineWidth',0.5, 'Linecolor', 'black')
     hold off
     xlim(100*[0,L]), ylim(100*[0,L])
     title('V_D aktiv')
     xlabel('x in cm'), ylabel('y in cm')
     set(gca,'FontSize',18)
     
 figure()
 set(gcf, 'Position', get(0, 'Screensize'));
 [gx,gy] = gradient(-(V1+V2+V3+V4),hx,hy);
 quiver(100*X(1:size(X,1)/25:end,1:size(X,2)/25:end),...
     100*Y(1:size(Y,1)/25:end,1:size(Y,2)/25:end),...
     gx(1:size(gx,1)/25:end,1:size(gx,2)/25:end),...
     gy(1:size(gy,1)/25:end,1:size(gy,2)/25:end),100);
 hold on, contour(X*100,Y*100,V1+V2+V3+V4,'LineWidth',0.5, 'Linecolor', 'black')
 hold off
 xlim(100*[0,L]), ylim(100*[0,L])
 title('Alle V_i aktiv')
 xlabel('x in cm'), ylabel('y in cm')
 set(gca,'FontSize',18)
 
 
 %% Save the plots
dir_out = 'Bsp_1_analytical_figures';
[~,~] = mkdir(dir_out);
h =  findobj('type','figure');
for k = 1 : length(h)
    f = h(k);
    name = fullfile(dir_out,sprintf(['Task',task_var,'_fig_%d'],f.Number));
    %savefig(h(k),name);
    saveas(h(k),name,'epsc')
    
end

 
 %% Functions
 function [V1,V2,V3,V4] = calc_V (x,y,N_sum,L,VA,VB,VC,VD)
 
 V1 = 0;
 V2 = 0;
 V3 = 0;
 V4 = 0;
 for k = 1 : N_sum
     
     n=2*k-1;
     
     V1 = V1 + 4*VA/(n*pi) * sinh(n*pi*(L-y)/L)/sinh(n*pi) * sin(n*pi*x/L);
     V2 = V2 + 4*VB/(n*pi) * sinh(n*pi*x/L)/sinh(n*pi) * sin(n*pi*y/L);
     V3 = V3 + 4*VC/(n*pi) * sinh(n*pi*y/L)/sinh(n*pi) * sin(n*pi*x/L);
     V4 = V4 + 4*VD/(n*pi) * sinh(n*pi*(L-x)/L)/sinh(n*pi) * sin(n*pi*y/L);
 end
 
 end
 

