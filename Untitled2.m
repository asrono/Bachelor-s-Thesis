clear all; close all;
power = linspace(-1,4,100);
m = 100;
theta  = linspace(0,2*pi,m);
x = 1.1*exp(1i*theta);
figure(1)
plot(real(x),imag(x),'.k'); hold on
axis([-1,1,-1,1]*2);
grid on
title('z')
for n = power
y = log(x);
figure(2)
plot(real(y),imag(y),'.k'); hold on
axis([-1,1,-1,1]*4);
grid on
title(strcat('f(z) = z^{',strcat(string(n)),'}'));
hold off
end