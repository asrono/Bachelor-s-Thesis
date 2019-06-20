clear all; close all;

figure()
x  = linspace(-1,1,1e3);

y = 1+0.5*cos(x*2*pi);

plot(x,y); hold on;
plot(x,sin(x)./(x).*y)
