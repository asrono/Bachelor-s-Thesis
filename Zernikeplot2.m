clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex');

a = 1e2;
r = linspace(0,1,a);
theta = linspace(0,2*pi,a);
n_max = 4;
figure(1);
ind = 1;
for i = 1:(n_max+1)
    n = i - 1;
    m_arr = -n:2:n;
    for j = 1:length(m_arr)
        m = m_arr(j);
        subplot(n_max+1,n_max+1,ind)
            Z = Zer(n,m,r,theta);
            X = r.*sin(theta)';
            Y = r.*cos(theta)';
            
            surf(X,Y,Z)
            colormap copper
            shading interp
            xlim([-1,1])
            ylim([-1,1])
            view(45,90)
            grid off
            axis off       
        ind = ind + 1;
    end
    ind = ind + n_max - length(m_arr)+1;
end

%%
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