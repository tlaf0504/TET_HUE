clear all
close all
clc

V1 = 1;
V2 = 1;
V3 = 1;
V4 = 1;


N1 = @(xi, eta) (1-xi).*(1-eta);
N4 = @(xi, eta) xi.*(1-eta);
N3 = @(xi, eta) xi.*eta;
N2 = @(xi, eta) (1-xi).*eta;

xi = linspace(0,1,100);
eta = linspace(0,1,100);
zeros_vec = zeros(length(eta), 1);
ones_vec = ones(length(eta), 1);
element_outlines_x = [0,1,1,0,0];
element_outlines_y = [0,0,1,1,0];
element_outlines_z = [0,0,0,0,0];

figure()
p1 = subplot(221);
plot3(element_outlines_x, element_outlines_y, element_outlines_z,'Color', 'Black')
hold on
plot3(xi,zeros_vec,N1(xi,0), 'Color', 'red')
plot3(zeros_vec,eta,N1(0,eta), 'Color', 'green')
plot3(xi,eta,N1(xi,eta), 'Color', 'blue')
grid minor
title('N_1(\xi, \eta)')
xlabel('\xi')
ylabel('\eta')


p2 = subplot(222);
plot3(element_outlines_x, element_outlines_y, element_outlines_z,'Color', 'Black')
hold on
plot3(zeros_vec,1-eta,N2(0,1-eta), 'Color', 'red')
plot3(xi, ones_vec, N2(xi,1), 'Color', 'green')
plot3(xi,1-eta,N2(xi,1-eta), 'Color', 'blue')
grid minor
title('N_2(\xi, \eta)')
xlabel('\xi')
ylabel('\eta')


p3 = subplot(223);
plot3(element_outlines_x, element_outlines_y, element_outlines_z,'Color', 'Black')
hold on
plot3(1-xi,ones_vec,N3(1-xi,1), 'Color', 'red')
plot3(ones_vec,1-eta,N3(1,1-eta), 'Color', 'green')
plot3(1-xi,1-eta,N3(1-xi,1-eta), 'Color', 'blue')
grid minor
title('N_3(\xi, \eta)')
xlabel('\xi')
ylabel('\eta')


p4 = subplot(224);
plot3(element_outlines_x, element_outlines_y, element_outlines_z,'Color', 'Black')
hold on
plot3(1-xi,zeros_vec,N4(1-xi,0), 'Color', 'red')
plot3(ones_vec,eta,N4(1,eta), 'Color', 'green')
plot3(1-xi,eta,N4(1-xi,eta), 'Color', 'blue')
grid minor
title('N_4(\xi, \eta)')
xlabel('\xi')
ylabel('\eta')

