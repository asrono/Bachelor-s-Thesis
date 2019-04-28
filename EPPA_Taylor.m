%% EPPA_Taylor
% This code makes a plots to demonstrate taylor expansions.

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 28-04-2019

% Detailed description:

%% Settings
clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex');

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\EPPA\Taylor/';

xi = -2.5; xf = -xi;

x_ini = linspace(xi,xf,25);
y = sin(x_ini);

x = linspace(xi*1.1,xf*1.1,1e3);

y1 = x;
y2_1 = x.^2;
y2_2 = x.^3;
y2_3 = -x.^3;
y2 = -x.^3/factorial(3);
y3_1 = x.^4;
y3 = x.^5/factorial(5);

figure(1)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
legend({'$f(x)$'}, 'interpreter','latex','location','southeast')
grid on
figure_name = 'EPPA_Taylor_1';
filetype    = '.png';
print(figure(1), '-dpng', strcat(folder,figure_name,filetype))

figure(2)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1,'r-');
legend({'$f(x)$','$x$'}, 'interpreter','latex','location','southeast')
figure_name = 'EPPA_Taylor_2';
filetype    = '.png';
print(figure(2), '-dpng', strcat(folder,figure_name,filetype))

figure(3)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2_1,'r-');
legend({'$f(x)$','$x+x^2$'}, 'interpreter','latex','location','southeast')
figure_name = 'EPPA_Taylor_3';
filetype    = '.png';
print(figure(3), '-dpng', strcat(folder,figure_name,filetype))

figure(4)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2_2,'r-');
legend({'$f(x)$','$x+x^3$'}, 'interpreter','latex','location','southeast')
figure_name = 'EPPA_Taylor_4';
filetype    = '.png';
print(figure(4), '-dpng', strcat(folder,figure_name,filetype))

figure(5)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2_3,'r-');
legend({'$f(x)$','$x-x^3$'}, 'interpreter','latex','location','southeast')
figure_name = 'EPPA_Taylor_5';
filetype    = '.png';
print(figure(5), '-dpng', strcat(folder,figure_name,filetype))

figure(6)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2,'r-');
legend({'$f(x)$','$x-\frac{x^3}{6}$'}, 'interpreter','latex','location','southeast')
figure_name = 'EPPA_Taylor_6';
filetype    = '.png';
print(figure(6), '-dpng', strcat(folder,figure_name,filetype))

figure(7)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2+y3_1,'r-');
legend({'$f(x)$','$x-\frac{x^3}{6}+x^4$'}, 'interpreter','latex','location','southeast')
figure_name = 'EPPA_Taylor_7';
filetype    = '.png';
print(figure(7), '-dpng', strcat(folder,figure_name,filetype))

figure(8)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2+y3,'r-');
legend({'$f(x)$','$x-\frac{x^3}{6}+\frac{x^5}{120}$'}, 'interpreter','latex','location','southeast')
figure_name = 'EPPA_Taylor_8';
filetype    = '.png';
print(figure(8), '-dpng', strcat(folder,figure_name,filetype))