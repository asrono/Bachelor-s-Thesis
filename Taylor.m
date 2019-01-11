clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex');

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
print('EPPA_Taylor_1','-dpng','-r300');

figure(2)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1,'r-');
legend({'$f(x)$','$x$'}, 'interpreter','latex','location','southeast')
print('EPPA_Taylor_2','-dpng','-r300');

figure(3)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2_1,'r-');
legend({'$f(x)$','$x+x^2$'}, 'interpreter','latex','location','southeast')
print('EPPA_Taylor_3','-dpng','-r300');

figure(4)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2_2,'r-');
legend({'$f(x)$','$x+x^3$'}, 'interpreter','latex','location','southeast')
print('EPPA_Taylor_4','-dpng','-r300');

figure(5)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2_3,'r-');
legend({'$f(x)$','$x-x^3$'}, 'interpreter','latex','location','southeast')
print('EPPA_Taylor_5','-dpng','-r300');

figure(6)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2,'r-');
legend({'$f(x)$','$x-\frac{x^3}{6}$'}, 'interpreter','latex','location','southeast')
print('EPPA_Taylor_6','-dpng','-r300');

figure(7)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2+y3_1,'r-');
legend({'$f(x)$','$x-\frac{x^3}{6}+x^4$'}, 'interpreter','latex','location','southeast')
print('EPPA_Taylor_7','-dpng','-r300');

figure(8)
plot(x_ini,y,'ko');hold on
xlim([xi,xf]*1.1);
ylim([-1,1]*2)
grid on
plot(x,y1+y2+y3,'r-');
legend({'$f(x)$','$x-\frac{x^3}{6}+\frac{x^5}{120}$'}, 'interpreter','latex','location','southeast')
print('EPPA_Taylor_8','-dpng','-r300');