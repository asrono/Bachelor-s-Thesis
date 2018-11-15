clear all; close all;
power = linspace(-1,-1,500);
m = 500;

x_1 = linspace(0,1,m)*1i;
x_2 = linspace(0,1,m);
x_3 = linspace(0,1,m)*1i+1;
x_4 = linspace(0,1,m)+1i;
figure(1)
plot(real(x_1),imag(x_1),'.g'); hold on
plot(real(x_2),imag(x_2),'.b');
plot(real(x_3),imag(x_3),'.r');
plot(real(x_4),imag(x_4),'.y');
axis([-1,1,-1,1]*2);
grid on
title('z')
for n = power
y_1 = x_1.^n;
y_2 = x_2.^n;
y_3 = x_3.^n;
y_4 = x_4.^n;
figure(2)
plot(real(y_1),imag(y_1),'.g'); hold on
plot(real(y_2),imag(y_2),'.b');
plot(real(y_3),imag(y_3),'.r');
plot(real(y_4),imag(y_4),'.y');
axis([-1,1,-1,1]*4);
grid on
title(strcat('f(z) = z^{',strcat(string(n)),'}'));
hold off
end