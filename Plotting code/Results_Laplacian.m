%% Results_Laplacian
% Checks the validity of the Laplacian expressions found

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

%% Creating plotting vectors
nplot = 5e3;
nplot2 = 5e1;
h = 1 / (nplot - 1);
rplot = linspace(0,1,nplot);
rplot2 = linspace(0.001,1,nplot2);

surface = coef2surf(a_vec,rplot);
surface_num = diff(surface,2)/h^2 + 1 ./ rplot(1:end-2) .* diff(surface(1:end-1))/h;
surface_sum = coef2surf_laplacian(a_vec,rplot2);
surface_ana = coef2surf_laplacian2(a_vec,rplot2);
surface_jan = coef2surf(Zernike_laplacian_Janssen_vec(a_vec),rplot2);
surface_jan_inv = coef2surf(Zernike_inverse_laplacian_Janssen_vec(Zernike_laplacian_Janssen_vec(a_vec)),rplot2);
%% Plotting
figure(1); title('Expression for Laplacian');
subplot(2,2,1); hold on;
title('Surface');
plot(rplot,surface,'-k');
plot(rplot2,surface_jan_inv,'ok');
xlabel('$r$');
ylabel('$h(r)$')
grid on
legend({'Original surface','Surface via Janssen'},...
        'Interpreter','latex',...
        'Location', 'northeast');
  
subplot(2,2,2); hold on;
title('Laplacian of surface');
plot(rplot(1:end-2),surface_num,'-k');
plot(rplot2,surface_sum,'ok');

grid on
ylim([-2,4])
xlim([0,1])
xlabel('$r$')
ylabel('$\bigtriangleup r$')
legend({'Numeric','Summation'},...
        'Interpreter','latex',...
        'Location', 'northeast');
    
subplot(2,2,3); hold on;
title('Surface');
plot(rplot(1:end-2),surface_num,'-k');
plot(rplot2,surface_ana,'ok');     
grid on
ylim([-2,4])
xlim([0,1])
xlabel('$r$')
ylabel('$\bigtriangleup r$')    
legend({'Numeric','Novel expression'},...
        'Interpreter','latex',...
        'Location', 'northeast');       

subplot(2,2,4); hold on;
title('Surface');
plot(rplot(1:end-2),surface_num,'-k');
plot(rplot2,surface_jan,'ok');
grid on
ylim([-2,4])
xlim([0,1])
xlabel('$r$')
ylabel('$\bigtriangleup r$')
legend({'Numeric','Janssen'},...
        'Interpreter','latex',...
        'Location', 'northeast');    

%% End
disp('Done!')
%% Functions
function a_vec_out = Zernike_inverse_laplacian_Janssen_vec(a_vec_in)
    a_vec_out = zeros(1,length(a_vec_in)+100);
    for i = 1: length(a_vec_in)
        if a_vec_in(i) == 0
           continue 
        end
        [n,m] = single2double_index(i-1);
        A = Zernike_inverse_laplacian_Janssen(n,m);
        A(:,3) = a_vec_in(i) * A(:,3);
        terms = length(A(:,3));
        for k = 1:terms
            j = double2single_index(A(k,1),A(k,2));
            a_vec_out(j+1) = a_vec_out(j+1) + A(k,3);
        end
    end
end

function A = Zernike_inverse_laplacian_Janssen(n,m)
    m_new = m*ones(1,3);
    n_new = (n-2):2:(n+2);   
    coef = [1/(4*n*(n+1)),...
            -1/(2*n*(n+2)),...
            1/(4*(n+1)*(n+2))];
    
    A(1,:) = n_new;
    A(2,:) = m_new;
    A(3,:) = coef;
    
    check = A(1,:) - abs(A(2,:));
    logis = check>0;
    A(3,:) = coef.*logis;
    
    A(isnan(A))=0;
    
    logis2 = ~any(A(3,:),1);
    A(:,logis2) = [];
    
    A = A';
end


function a_vec_out = Zernike_laplacian_Janssen_vec(a_vec_in)
    a_vec_out = zeros(1,length(a_vec_in));
    for i = 1: length(a_vec_in)
        if a_vec_in(i) == 0
           continue 
        end
        [n,m] = single2double_index(i-1);
        if n == m
            continue
        end
        A = Zernike_laplacian_Janssen(n,m);
        A(:,3) = a_vec_in(i) * A(:,3);
        terms = length(A(:,3));
        for k = 1:terms
            j = double2single_index(A(k,1),A(k,2));
            a_vec_out(j+1) = a_vec_out(j+1) + A(k,3);
        end
    end
end

function A = Zernike_laplacian_Janssen(n,m)
    s_list = abs(m):2:(n-2);
    terms = length(s_list);
    m_new = m*ones(1,terms);
    n_new = zeros(1,terms);
    coef  = zeros(1,terms);

    index = 1;
    for s = s_list
        n_new(index) = s;
        coef(index)  = (s+1)*(n+s+2)*(n-s);
        
        index = index + 1;
    end
    A(1,:) = n_new;
    A(2,:) = m_new;
    A(3,:) = coef;
    
    A = A';
end

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
