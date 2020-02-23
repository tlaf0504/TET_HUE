clear all
close all
clc

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

 
 figure()
 p1 = subplot(221);
 surf(X,Y,V1,'EdgeColor','none');
 colorbar
 colormap jet
 title('V_1')
 
 p2 = subplot(222);
 surf(X,Y,V2,'EdgeColor','none');
 colorbar
 colormap jet
 title('V_2')
 
 p3 = subplot(223);
 surf(X,Y,V3,'EdgeColor','none');
 colorbar
 colormap jet
 title('V_3')
 
 p4 = subplot(224);
 surf(X,Y,V4,'EdgeColor','none');
 colorbar
 colormap jet
 title('V_4')
 
 figure()
 surf(X,Y,V1+V2+V3+V4,'EdgeColor','none');
 colorbar
 colormap jet
 title('V_{ges}')
 
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
