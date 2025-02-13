%% Laplacian
% Creates ray tracer for Zernike purposes and laplacian magic window

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 21-05-2019

% Detailed description:

%% Settings
close all; clear all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex'); 

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\Ray Tracer Zernike\';

%% Initialise systema

% At x = xi, we place a lens with known zernike polynomials.
% At x = xf, we place a screen and observe intensity.

xi = 0; % arbitrary units
xf = 1/2; % arbitrary units

lens = [-1,1]; % position of bottom and top of lens
lens_middle = lens(1) + ( lens(2)-lens(1) );


%% Settings


%% Input
% load coefficients for a lens
% loadfile = 'a_vec_Berry.mat';
loadfile = 'a_vec_Berry.mat';

load(strcat('C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Data\',loadfile),'a_vec_N');
a_vec = a_vec_N;    % Rename variable for this script
a_vec(1) = xi;       % Set height of lens to optical start
a_vec = a_vec/2;
% a_vec = zeros(1,25);
% a_vec(9) = 0.0001;
% a_vec(13) = 0.02;
% a_vec(25) = 0.01;
clear a_vec_N

%% Plotting
figure(1); title('Correctness of the Laplacian'); hold on;
nplot = 1e3; h = 1 / (nplot - 1);
rplot = linspace(0,1,nplot);

surface = coef2surf(a_vec,rplot);
plot(rplot,200*surface,'-k');

surface2 = coef2surf_deriv(a_vec,rplot);
plot(rplot,surface2,'-r');

surface3 = coef2surf_laplacian(a_vec,rplot);
plot(rplot,surface3,'-r','LineWidth',5);

surface4 = coef2surf_laplacian2(a_vec,rplot);
plot(rplot,surface4,'-g','LineWidth',2);

surface5 = diff(surface,2)/h^2 + 1 ./ rplot(1:end-2) .* diff(surface(1:end-1))/h;
plot(rplot(1:end-2),surface5,':k','LineWidth',2);
ylim([-2,4])
legend({'Surface','Derivative','Laplacian (analytical different basis)','Laplacian (analytical same basis)','Laplacian Numerical'},...
        'Interpreter','latex',...
        'Location', 'northeast');
    
%% Find Intensity
z = linspace(0,0.5,1e3);
n = 1.6;
I = 1 - (n - 1)*z'.*coef2surf_laplacian2(a_vec,rplot);

I2 = 1./ (1 + (n - 1)*z'.*coef2surf_laplacian(a_vec,rplot));
figure(2);
plot(rplot,I(end,:),'g'); hold on;
plot(rplot,I2(end,:),'r');
plot(rplot,I2(end,:)-I(end,:),'k');
legend({'linear approx','non linear approx','diff'},...
        'Interpreter','latex',...
        'Location', 'northeast');

figure(3);
surf(I);
shading interp

figure(4);
surf(I2);
shading interp

figure(1)
%% End
disp('Done!')
%% Functions
function Z = Zernike_surface(a_vec,r,theta)
    Z = zeros(numel(r),numel(theta))';
    i = 0;
    for a = a_vec
        [n,m] = single2double_index(i);
        Z = Z + a*Zer(n,m,r,theta);
        i = i + 1;
    end
end

function rad = coef2surf(a_vec,r_plot)
    for j = 0:(length(a_vec)-1)
        [n,m] = single2double_index(j);
        rad_mat(j+1,:) = a_vec(j+1)*R(n,m,r_plot);
    end
    rad = sum(rad_mat,1);
end

function rad = coef2surf_deriv(a_vec,r_plot)
    for j = 0:(length(a_vec)-1)
        [n,m] = single2double_index(j);
        rad_mat(j+1,:) = a_vec(j+1)*R_deriv(n,m,r_plot);
    end
    rad = sum(rad_mat,1);
end

function rad = coef2surf_deriv_num(a_vec,r_plot)
    for j = 0:(length(a_vec)-1)
        [n,m] = single2double_index(j);
        rad_mat(j+1,:) = a_vec(j+1)*R_deriv_num(n,m,r_plot);
    end
    rad = sum(rad_mat,1);
end

function rad = coef2surf_deriv2_num(a_vec,r_plot)
    for j = 0:(length(a_vec)-1)
        [n,m] = single2double_index(j);
        rad_mat(j+1,:) = a_vec(j+1)*R_deriv2_num(n,m,r_plot);
    end
    rad = sum(rad_mat,1);
end

function rad = coef2surf_laplacian(a_vec,rplot)
    rad = coef2surf_deriv2_num(a_vec,rplot) + 1./rplot .* coef2surf_deriv_num(a_vec,rplot);
end

function rad = coef2surf_laplacian2(a_vec,r_plot)
    for j = 0:(length(a_vec)-1)
        [n,m] = single2double_index(j);
        rad_mat(j+1,:) = a_vec(j+1)*R_laplacian(n,m,r_plot);
    end
    rad = sum(rad_mat,1);
end
function [surface, r, r_max] = extract_points(lens)
    surface(1,:) = lens(1,:,:); % get scalar values for lens surface
    r_max = max(lens(2,:,:));   % find maximum radius element
    r(1,:) = lens(2,:,:)/r_max; % get normalised radius points
end

function [r_interpolated, surface_interpolated] = ...
    interpol(surface, r,r_int_min,r_int_max,n_interpolation)

    r_interpolated = linspace(r_int_min,r_int_max,n_interpolation);
    surface_interpolated = interp1(r,surface,r_interpolated,'spline');    
end

function in = extend_to_neg(in)
    in2 = in;
    in2(2,:,:) = -in(2,:,:);
    in(:,:,(end+1):(2*length(in))) = in2;
end

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

function radial = R_deriv(n,m,r)
    radial = (r.^2.*(n+2)+m)./(r.*(1-r.^2)).*R(n,m,r) - (n+m+2)./(1-r.^2).*R(n+1,m+1,r);
end

function radial = R_laplacian(n,m,r)
%     c0 = 1./(r.^2.*(1-r.^2).^2);
% %     c1 = r.^6*(n+8) + r.^4*(n+3*m+2+(n+2)^2) + r.^2*(n+3*m+2*n*m+2) + m*(m+1);
%     c1 = r.^4*(n+2)*(n+1) + r.^2*((n+2)*(2+2*m)+3*m) + m^2;
%     c2 = -(n+m+2)*(-r.^3 + r.^2*(2*n+7) + r + 2*m+1);
%     c3 = (n+m+2)*r.^2*(n+m+4);
%     radial = c0 .* (c1.*R(n,m,r) + c2.*R(n+1,m+1,r) + c3.*R(n+2,m+2,r));

%     radial = (m*(3*r.^2-1)+(n+2)*(r.^2+1).*r.^2)./(r.^2.*(1-r.^2).^2).*R(n,m,r)+...
%         (r.^2*(n+2)+m)./(r.*(1-r.^2)).*((r.^2*(n+2)+m)./(r.*(1-r.^2)).*R(n,m,r)-(n+m+2)./(1-r.^2).*R(n+1,m+1,r))-...
%         2*r*(n+m+2)./(1-r.^2).^2.*R(n+1,m+1,r) -...
%         (n+m+2)./(1-r.^2).*((r.^2*(n+3)+m+1)./(r.*(1-r.^2)).*R(n+1,m+1,r)-(n+m+4)./(1-r.^2).*R(n+2,m+2,r))+...
%         ((r.^2*(n+2)+m)./(r.^2.*(1-r.^2)).*R(n,m,r)-(n+m+2)./(r.*(1-r.^2)).*R(n+1,m+1,r));

%     c0 = 1./(r.^2.*(1-r.^2).^2);
%     radial = c0 .* (...
%         R(n  ,m  ,r).*(m*(3*r.^2-1)+(n+2)*(1+r.^2).*r.^2+(r.^2*(n+2)+m).^2)-...
%         R(n+1,m+1,r).*(r.*(r.^2*(n+2)+m)*(n+m+2)+2*r.^3*(n+m+2)+r.*(n+m+2).*(r.^2*(n+3)+m+1))+...
%         R(n+2,m+2,r).*(r.^2*(n+m+2)*(n+m+4)) )+...
%         ((r.^2*(n+2)+m)./(r.^2.*(1-r.^2)).*R(n,m,r)-(n+m+2)./(r.*(1-r.^2)).*R(n+1,m+1,r));

%     c0 = 1./(r.^2.*(1-r.^2).^2);
%     radial = c0 .* (...
%         R(n  ,m  ,r).*(r.^4*(n+2)*(n+3)+r.^2*(m*(2*n+7)+n+2)+m*(m-1))-...
%         R(n+1,m+1,r).*((n+m+2)*r.*(r.^2*(2*n+7)+2*m+1))+...
%         R(n+2,m+2,r).*(r.^2*(n+m+2)*(n+m+4)) )+...
%         ((r.^2*(n+2)+m)./(r.^2.*(1-r.^2)).*R(n,m,r)-(n+m+2)./(r.*(1-r.^2)).*R(n+1,m+1,r));
   
%     c0 = 1./(r.^2.*(1-r.^2).^2);
%     radial = c0 .* (...
%         R(n  ,m  ,r).*(r.^4*(n+2)*(n+3)+r.^2*(m*(2*n+7)+n+2)+m*(m-1)+(1-r.^2).*(r.^2*(n+2)+m))-...
%         R(n+1,m+1,r).*((n+m+2)*r.*(r.^2*(2*n+7)+2*m+1+(1-r.^2)))+...
%         R(n+2,m+2,r).*(r.^2*(n+m+2)*(n+m+4)) );

%     Last Stable Version
    c0 = 1./(r.^2.*(1-r.^2).^2);
    radial = c0 .* (...
        R(n  ,m  ,r).*(r.^4*(n+2)^2+2*r.^2*(m*(n+3)+n+2)+m^2)-...
        R(n+1,m+1,r).*(n+m+2).*r.*(r.^2*(2*n+6)+2*m+2)+...
        R(n+2,m+2,r).*(r.^2*(n+m+2)*(n+m+4)) );


%     c0 = 1./(r.^2.*(1-r.^2).^2*(n+2));
%     radial = c0 .* (...
%         R(n  ,m  ,r).*(r.^4*(n^3+6*n^2+12*n+8)+r.^2*(n*m^2+3*m^2+2*m*n^2+12*m*n+18*m-n^3-3*n^2+2*n+8)+m^3+m^2*n+5*m^2-m*n^2-2*m*n+2*m-n^2-2*n)-...
%         R(n+2,m+2,r)*(n+m+2)*(n+m+4).*(m+r.^2+1));    

%     c0 = 1./(r.^2.*(1-r.^2).^2);
%     radial = c0 .* (...
%         R(n  ,m  ,r).*(r.^4*(n+2)^2+2*r.^2*(m*(n+3)+n+2)+m^2-(n-m)/(2*(n+2))*(n+m+2)*(r.^2*(2*n+6)+2*m+1))+...
%         R(n+2,m+2,r).*(r.^2*(n+m+2)*(n+m+4) - (n+m+4)/(2*(n+2))*(n+m+2)*(r.^2*(2*n+6)+2*m+1))); 
end

function radial = R_deriv_num(n,m,r)
    m = abs(m);
    top = (n-m)/2;
    bot = (n+m)/2;
    for s = 0 : top
        x(s+1,:) =  (-1)^s*factorial(n-s)/(factorial(s)*factorial(bot-s)*factorial(top-s))*(n-2*s)*r.^(n-2*s-1);
    end
    if top ~= 0
    radial = sum(x);
    else
    radial = x;
    end
end

function radial = R_deriv2_num(n,m,r)
    m = abs(m);
    top = (n-m)/2;
    bot = (n+m)/2;
    for s = 0 : top
        x(s+1,:) =  (-1)^s*factorial(n-s)/(factorial(s)*factorial(bot-s)*factorial(top-s))*(n-2*s)*(n-2*s-1)*r.^(n-2*s-2);
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

function a = zernikecoef(n,m,r,f)
    norm = 1/(2*(n+1));
    int = trapz(r,R(n,m,r).*f.*r);
    a = int/norm;
end
function plot_point(P,c)
    plot(P(1),P(2),c)
end

function plot_points(P_vec,c)
    hold on;
    for P = P_vec       
        plot_point(P,c);
    end
    hold off;
end

function j = double2single_index(n,m)
    j = ( n * (n+2)+ m) / 2;
end

function [n,m] = single2double_index(j)
    n = ceil( (-3 + (9 + 8*j)^(1/2)) / 2 );
    m = 2*j - n * (n+2);
end

function v_r = find_reflected_ray(v_i,n_i,n_r,normal)
    if dot(v_i,normal) < 0
        normal = -normal;
        c = -1;
    else
        c = 1;
    end
    n = n_i/n_r;
    v_r = n*v_i+(n*dot(normal,v_i)-sqrt(1-n^2*(1-dot(v_i,normal)))) * normal*c;
    v_r = v_r./norm(v_r);
end
