clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',18);
set(0,'defaultAxesTickLabelInterpreter','latex');

a = 1e3;
r = linspace(0,1,a);
theta = linspace(0,2*pi,a);

n=4;
m=-2;

Z = Zer(n,m,r,theta);
X = r.*sin(theta)';
Y = r.*cos(theta)';

title_text = strcat('$n=',string(n),', m=',string(m),'$');


surf(X,Y,Z)
colormap jet
shading interp
xlim([-1,1])
ylim([-1,1])
xlabel('$x$');ylabel('$y$');zlabel('$z$');
view(45,45)
title(title_text)
grid on
axis tight


filename = strcat('Single_Zernike_',string(n),string(m));
print(filename,'-dpng','-r300');

function radial = R(n,m,r)
    m = abs(m);
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

function zernike = Zer(n,m,r,theta)
    if m >= 0
        zernike = R(n,m,r).*cos(m*theta)';
    else
        zernike = R(n,m,r).*sin(m*theta)';
    end
end