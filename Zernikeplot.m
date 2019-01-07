clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex');  

r = 0:1e-3:1;
%%
m = 0;
n_lst = [2,4,6,8];
set = {'k-','k:','k-.','k--'};
figure('position', [1000, 700, 1000, 500])
for i = 1:length(n_lst)
    n = n_lst(i);
    y = R(n,m,r);
    plot(r,y,set{i}); hold on
end
hold off;
title('Zernike even radial function')
xlabel('$r$'); ylabel('$R_n^0(r)$');
xlim([0,1])
ylim([-1,1])
grid on
legend({'$n=2$','$n=4$','$n=6$','$n=8$'},'interpreter','latex','location','southeast')
pbaspect([2,1,1])
print('Zernike_radial_even','-dpng','-r300');
%%
m = 1;
n_lst = [1,3,5,7];
set = {'k-','k:','k-.','k--'};
figure('position', [1000, 100, 1000, 500])
for i = 1:length(n_lst)
    n = n_lst(i);
    y = R(n,m,r);
    plot(r,y,set{i}); hold on
end
hold off;
title('Zernike odd radial function')
xlabel('$r$'); ylabel('$R_n^1(r)$');
xlim([0,1])
ylim([-1,1])
grid on
legend({'$n=1$','$n=3$','$n=5$','$n=7$'},'interpreter','latex','location','southeast')
pbaspect([2,1,1])
print('Zernike_radial_odd','-dpng','-r300');
%%
n = 7;
m_lst = [1,3,5,7];
set = {'k-','k:','k-.','k--'};
figure('position', [0, 700, 1000, 500])
for i = 1:length(m_lst)
    m = m_lst(i);
    y = R(n,m,r);
    plot(r,y,set{i}); hold on
end
hold off;
title('Zernike radial function')
xlabel('$r$'); ylabel('$R_7^m(r)$');
xlim([0,1])
ylim([-1,1])
grid on
legend({'$m=1$','$m=3$','$m=5$','$m=7$'},'interpreter','latex','location','southeast')
pbaspect([2,1,1])
print('Zernike_radial_m','-dpng','-r300');

%%
function radial = R(n,m,r)
    top = (n-m)/2;
    bot = (n+m)/2;
    for s = 0 : top
        x(s+1,:) =  (-1)^s*factorial(n-s)/(factorial(s)*factorial(bot-s)*factorial(top-s))*r.^(n-2*s);
    end
    if top ~= 0
    radial = sum(x);
    else
    radial = x;
    end
end