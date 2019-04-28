clear all; close all;set(0,'defaulttextinterpreter','latex'); set(0,'defaultaxesfontsize',14);
% rgb = 1 is color, rgb = 0 is black and white
rgb = 1;
% define lower and upper bound for power 
lower = -10;
upper = 13;
power = linspace(lower,upper,1000);
% define number of points per line segment
m = 1500;
% define linewidth/markersize
o = 5;
% define line segments
x_1 = linspace(-1,1,m)*1i-1;
x_2 = linspace(-1,1,m)-1i;
x_3 = linspace(-1,1,m)*1i+1;
x_4 = linspace(-1,1,m)+1i;
% plot pre-image
figure(1);
set(groot,'defaultLineLineWidth',2)
movegui('west');
if rgb
plot(real(x_1),imag(x_1),'.b','MarkerSize',o); hold on
plot(real(x_2),imag(x_2),'.','Color',[0.996, 0.6400, 0.00],'MarkerSize',o);
plot(real(x_3),imag(x_3),'.g','MarkerSize',o);
plot(real(x_4),imag(x_4),'.r','MarkerSize',o); hold off
else
plot(real(x_1),imag(x_1),'.k','MarkerSize',o); hold on
plot(real(x_2),imag(x_2),'.k','MarkerSize',o);
plot(real(x_3),imag(x_3),'.k','MarkerSize',o);
plot(real(x_4),imag(x_4),'.k','MarkerSize',o); hold off
end
axis([-1,1,-1,1]*2);
set(gca,'ticklabelinterpreter','latex');
xlabel('Re');
ylabel('Im');
title('$z$');
% plot image
figure(2);
set(groot,'defaultLineLineWidth',2)
movegui('center')
for n = power
y_1 = x_1.^n;
y_2 = x_2.^n;
y_3 = x_3.^n;
y_4 = x_4.^n;
figure(2);
if rgb
plot(real(y_1),imag(y_1),'.b','MarkerSize',o); hold on
plot(real(y_2),imag(y_2),'.','Color',[0.996, 0.6400, 0.00],'MarkerSize',o);
plot(real(y_3),imag(y_3),'.g','MarkerSize',o);
plot(real(y_4),imag(y_4),'.r','MarkerSize',o); hold off
else
plot(real(y_1),imag(y_1),'.k','MarkerSize',o); hold on
plot(real(y_2),imag(y_2),'.k','MarkerSize',o);
plot(real(y_3),imag(y_3),'.k','MarkerSize',o);
plot(real(y_4),imag(y_4),'.k','MarkerSize',o); hold off
end
% change axis
if n < 0
    r = 1.1;
else
    r = 4;
end
axis([-1,1,-1,1]*r);
set(gca,'ticklabelinterpreter','latex');
xlabel('Re');
ylabel('Im');
% change title
title(strcat('$f(z) = z^{',strcat(string(n)),'}$'));
end