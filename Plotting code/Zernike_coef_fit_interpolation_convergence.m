%% Zernike coef fit
% This code fits zernike polynomials through data and outputs the coefs.
% Also plots the convergence of the coefficients.

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 28-04-2019

% Detailed description:

%% Settings
clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex');

%% Load data
savefile = 'LensExample.mat';
load(strcat('C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Data/',savefile));
l_OA = ER(1,1,3); % length of optical axis

%% Input handling
% update length of N
l_N = length(N_lens);

%% Fitting
% create empty vectors
surface = zeros(1,l_N);
r       = zeros(1,l_N);

% The zernike polynomials are defined on a disc with radius 1. Therefore we
% search for a represenation such that z = surface(r) such that |r| <= 1.
[surface_N, r_N, r_max_N] = extract_points(N_lens);
[surface_X, r_X, r_max_X] = extract_points(X_lens);

% plot raw datapoints
figure(2); movegui('northeast');
plot(r_N*r_max_N,surface_N,'ko')
hold on
plot(r_X*r_max_X,surface_X,'ko')

% set amount of interpolated points
n_list = 10.^linspace(log10(l_N),3,50);

% loop over n_list
for i = 1:length(n_list)
    n_interpolation = n_list(i);
    
    % get interpolated points
    [r_interpolated_N, surface_interpolated_N] = ...
        interpol(surface_N, r_N,0,1,n_interpolation);

    [r_interpolated_X, surface_interpolated_X] = ...
        interpol(surface_X, r_X,0,1,n_interpolation);

    % find coefficients
    m = 0;
    for n = 0:2:10
        a_vec_N(i,double2single_index(n,m)+1) = zernikecoef(n,m,r_interpolated_N,surface_interpolated_N); %#ok<SAGROW>
        a_vec_X(i,double2single_index(n,m)+1) = zernikecoef(n,m,r_interpolated_X,surface_interpolated_X); %#ok<SAGROW>
    end
    
    figure(1);
    semilogx(n_interpolation,a_vec_N(i,double2single_index(2,m)+1),'ok');
    hold on;
    semilogx(n_interpolation,a_vec_N(i,double2single_index(4,m)+1),'*k');
    semilogx(n_interpolation,a_vec_N(i,double2single_index(6,m)+1),'^k');
    semilogx(n_interpolation,a_vec_N(i,double2single_index(8,m)+1),'xk');
end
grid minor;
xlabel('$n_{\textrm{interpolation}}$');
ylabel('$a_{nm}$');
ylim([0,1.5]);
xlim([l_N,1e3]);
title('Convergence of coefficients')
legend({'$a_{20}$','$a_{40}$','$a_{60}$','$a_{80}$'},...
        'Interpreter','latex',...
        'Location', 'northeast');
movegui('east')
a_vec_N = a_vec_N(end,:);
a_vec_X = a_vec_X(end,:);

% find surface from zernike coef
r_plot = 0:0.01:1;
surface_zernike_N = coef2surf(a_vec_N,r_plot);
surface_zernike_X = coef2surf(a_vec_X,r_plot);

% plot zernike fit
figure(2);
plot(r_plot*r_max_N,surface_zernike_N,'-b')
plot(r_plot*r_max_X,surface_zernike_X,'-b')

% plot mesh
figure(3);
n_r = 1e1; % number mesh points radial direction
n_theta = 5e1; % number mesh points theta direction

r = linspace(0,1,n_r);
theta = linspace(0,2*pi,n_theta);

Z_N = a_vec_N(double2single_index(0,0)+1)*Zer(0,0,r,theta)+a_vec_N(double2single_index(2,0)+1)*Zer(2,0,r,theta)+a_vec_N(double2single_index(4,0)+1)*Zer(4,0,r,theta);
X_N = r_max_N*r.*sin(theta)';
Y_N = r_max_N*r.*cos(theta)';

Z_X = a_vec_X(double2single_index(0,0)+1)*Zer(0,0,r,theta)+a_vec_X(double2single_index(2,0)+1)*Zer(2,0,r,theta)+a_vec_X(double2single_index(4,0)+1)*Zer(4,0,r,theta);
X_X = r_max_X*r.*sin(theta)';
Y_X = r_max_X*r.*cos(theta)';

X = [X_N, fliplr(X_X)];
Y = [Y_N, fliplr(Y_X)];
Z = [Z_N, fliplr(Z_X)];
xlabel('$x$');ylabel('$y$');zlabel('$z$');
title('Mesh');
mesh(X,Y,Z);
shading interp
zlim([18,22]);
%%
function rad = coef2surf(a_vec,r_plot)
    for j = 0:(length(a_vec)-1)
        [n,m] = single2double_index(j);
        rad_mat(j+1,:) = a_vec(j+1)*R(n,m,r_plot);
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

function zernike = Zer(n,m,r,theta)
    if m >= 0
        zernike = R(n,m,r).*cos(m*theta)';
    else
        zernike = R(n,m,r).*sin(m*theta)';
    end
end

function a = zernikecoef(n,m,r,f)
    N_nm = 2*(1) / (1 + eq(m,0) );
    norm = trapz(r,r.*R(n,m,r).*R(n,m,r)*N_nm);
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