clear all; close all;
power = linspace(-5,5,200);

reals = linspace(0,1,50);
imags = transpose(linspace(0,1,50)*1i);

grid_points = reals + imags;

figure(1);
for n = power
plot(grid_points.^n); hold on;
plot( transpose(grid_points).^n ); hold off;
axis([-1,1,-1,1]*2);
grid on
title(strcat('f(z) = z*^{',strcat(string(n)),'}'));
drawnow;
end