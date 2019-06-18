%% Results Zernike Fit
% Plots figures to discuss goodness of fit

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 18-06-2019

% Detailed description:

%% Settings
close all; clear all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex'); 

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\Results\Zernike Fit\';

%% Set up variables
rplot = linspace(0,1,1e4);
rplot2 = linspace(0,1,20);

mu = 0;
sigma = 0.3;
y1 = normpdf(rplot,mu,sigma)*sqrt(2*pi*sigma^2);
y2 = 0.5+0.5*cos(rplot*2*pi);
y3 = 0.25+0.25*cos(rplot*4*pi)+normpdf(rplot,mu,sigma)*sqrt(1/2*pi*sigma^2);
y4 = 3/4*ones(1,numel(rplot)) - 1/2*[zeros(1,numel(rplot)/2) ones(1,numel(rplot)/2)];
y5 = 0.10+0.5*abs(cos(rplot*pi));
y6 = 0.25+0.5*sin(pi*rplot);

% find coefficients
m = 0;
for n = 0:2:10
    a_vec(double2single_index(n,m)+1) = zernikecoef(n,m,rplot,y1);
end
y1_fit1 = coef2surf(a_vec(1),rplot2);
y1_fit2 = coef2surf(a_vec(1:(double2single_index(2,0)+1)),rplot2);
y1_fit3 = coef2surf(a_vec(1:(double2single_index(6,0)+1)),rplot2);
y1_fit4 = coef2surf(a_vec(1:(double2single_index(10,0)+1)),rplot2);

% find coefficients
for n = 0:2:10
    a_vec2(double2single_index(n,m)+1) = zernikecoef(n,m,rplot,y2);
end
y2_fit1 = coef2surf(a_vec2(1),rplot2);
y2_fit2 = coef2surf(a_vec2(1:(double2single_index(2,0)+1)),rplot2);
y2_fit3 = coef2surf(a_vec2(1:(double2single_index(6,0)+1)),rplot2);
y2_fit4 = coef2surf(a_vec2(1:(double2single_index(10,0)+1)),rplot2);

% find coefficients
for n = 0:2:16
    a_vec3(double2single_index(n,m)+1) = zernikecoef(n,m,rplot,y3);
end
y3_fit1 = coef2surf(a_vec3(1),rplot2);
y3_fit2 = coef2surf(a_vec3(1:(double2single_index(10,0)+1)),rplot2);
y3_fit3 = coef2surf(a_vec3(1:(double2single_index(12,0)+1)),rplot2);
y3_fit4 = coef2surf(a_vec3(1:(double2single_index(16,0)+1)),rplot2);

% find coefficients
for n = 0:2:30
    a_vec4(double2single_index(n,m)+1) = zernikecoef(n,m,rplot,y4);
end
y4_fit1 = coef2surf(a_vec4(1),rplot);
y4_fit2 = coef2surf(a_vec4(1:(double2single_index(10,0)+1)),rplot);
y4_fit3 = coef2surf(a_vec4(1:(double2single_index(20,0)+1)),rplot);
y4_fit4 = coef2surf(a_vec4(1:(double2single_index(30,0)+1)),rplot);

% find coefficients
for n = 0:2:30
    a_vec5(double2single_index(n,m)+1) = zernikecoef(n,m,rplot,y5);
end
y5_fit1 = coef2surf(a_vec5(1),rplot);
y5_fit2 = coef2surf(a_vec5(1:(double2single_index(10,0)+1)),rplot);
y5_fit3 = coef2surf(a_vec5(1:(double2single_index(20,0)+1)),rplot);
y5_fit4 = coef2surf(a_vec5(1:(double2single_index(30,0)+1)),rplot);

% find coefficients
for n = 0:2:30
    a_vec6(double2single_index(n,m)+1) = zernikecoef(n,m,rplot,y6);
end
y6_fit1 = coef2surf(a_vec6(1),rplot);
y6_fit2 = coef2surf(a_vec6(1:(double2single_index(10,0)+1)),rplot);
y6_fit3 = coef2surf(a_vec6(1:(double2single_index(20,0)+1)),rplot);
y6_fit4 = coef2surf(a_vec6(1:(double2single_index(30,0)+1)),rplot);
%% Plotting
figure(1);
subplot(3,2,1);title('Surface 1: Gaussian'); hold on;
plot(rplot,y1,'-k');
plot(rplot2,y1_fit1,'xk')
plot(rplot2,y1_fit2,'ok')
plot(rplot2,y1_fit3,'*k')
plot(rplot2,y1_fit4,'sk')
legend({'Test function','$n = 0$','$n = 2$','$n = 6$','$n = 10$'},...
        'Interpreter','latex',...
        'Location', 'northeast');
grid on
xlabel('$r$')
ylabel('$h(r)$')
ylim([0,1])
% axis equal

subplot(3,2,2); title('Surface 2: cosine'); hold on;
plot(rplot,y2,'-k');
plot(rplot2,y2_fit1,'xk')
plot(rplot2,y2_fit2,'ok')
plot(rplot2,y2_fit3,'*k')
plot(rplot2,y2_fit4,'sk')
legend({'Test function','$n = 0$','$n = 2$','$n = 6$','$n = 10$'},...
        'Interpreter','latex',...
        'Location', 'north');
grid on
xlabel('$r$')
ylabel('$h(r)$')
ylim([0,1])



subplot(3,2,3); title('Surface 3: Gaussian + cosine'); hold on;
plot(rplot,y3,'-k');
plot(rplot2,y3_fit1,'xk')
plot(rplot2,y3_fit2,'ok')
plot(rplot2,y3_fit3,'*k')
plot(rplot2,y3_fit4,'sk')
legend({'Test function','$n = 0$','$n = 10$','$n = 12$','$n = 16$'},...
        'Interpreter','latex',...
        'Location', 'northeast');
grid on
xlabel('$r$')
ylabel('$h(r)$')
ylim([0,1])

subplot(3,2,4); title('Surface 4: discontinuity'); hold on;
plot(rplot,y4,'-k');
plot(rplot,y4_fit1,'--k')
plot(rplot,y4_fit2,'-.k')
plot(rplot,y4_fit3,':k')
plot(rplot,y4_fit4,'--k')
legend({'Test function','$n = 0$','$n = 10$','$n = 20$','$n = 30$'},...
        'Interpreter','latex',...
        'Location', 'northeast');
grid on
xlabel('$r$')
ylabel('$h(r)$')
ylim([0,1])

subplot(3,2,5); title('Surface 5: not continuously differentiable 1'); hold on;
plot(rplot,y5,'-k');
plot(rplot,y5_fit1,'--k')
plot(rplot,y5_fit2,'-.k')
plot(rplot,y5_fit3,':k')
plot(rplot,y5_fit4,'--k')
legend({'Test function','$n = 0$','$n = 10$','$n = 20$','$n = 30$'},...
        'Interpreter','latex',...
        'Location', 'north');
grid on
xlabel('$r$')
ylabel('$h(r)$')
ylim([0,1])

subplot(3,2,6); title('Surface 6: not continuously differentiable 2'); hold on;
plot(rplot,y6,'-k');
plot(rplot,y6_fit1,'--k')
plot(rplot,y6_fit2,'-.k')
plot(rplot,y6_fit3,':k')
plot(rplot,y6_fit4,'--k')
legend({'Test function','$n = 0$','$n = 10$','$n = 20$','$n = 30$'},...
        'Interpreter','latex',...
        'Location', 'south');
grid on
xlabel('$r$')
ylabel('$h(r)$')
ylim([0,1])
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
