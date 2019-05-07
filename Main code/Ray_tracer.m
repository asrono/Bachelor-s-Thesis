%% Ray tracer
% Creates ray tracer for Zernike purposes

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 07-05-2019

% Detailed description:

%% Settings
close all; clear all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex'); 

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\SMS Method/';

%% Initialise system

% At x = xi, we place a lens with known zernike polynomials.
% At x = xf, we place a screen and observe intensity.

xi = 0; % arbitrary units
xf = 10; % arbitrary units

lens = [-1,1]; % position of bottom and top of lens
lens_middle = ( lens(2)-lens(1) ) / 2;

%% Input
% load coefficients for a lens
loadfile = 'a_vec.mat';
load(strcat('C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Data\',loadfile),'a_vec_N');
a_vec = a_vec_N;    % Rename variable for this script
a_vec(1) = xi;       % Set height of lens to optical start
clear a_vec_N

%% Ray tracer
number_of_rays = 1e2;

% for plotting
y_plot = linspace(lens(1),lens(2),number_of_rays);
x_lens = coef2surf(a_vec,y_plot);

rays_start = linspace(lens(2),lens(1),number_of_rays)';

dy = zeros(number_of_rays,1);
dy(1) =   x_lens(2)     - x_lens(1);
dy(end) = x_lens(end-1) - x_lens(end);
for i = 2: (number_of_rays-1)
   dy(i) = ( x_lens(i+1)-x_lens(i-1) );
end
dx = ones(number_of_rays,1) * ( lens(1) - lens(2) ) / (number_of_rays - 1);
dx(1) = dx(1)/2;
dx(end) = dx(end)/2;
normal = [-dx,dy];
norm = sqrt( normal(:,1).^2 + normal(:,2).^2 );
normal = normal ./ norm;
%% plotting figure 1
figure(1); title('Optical system'); hold on; movegui('center')

% Plot optical axis
plot([xi,xf],lens_middle*ones(1,2),'-k')


% % Plot lens % placeholder
% plot(xi*ones(1,2),lens,':k')

% Plot lens (real lens)
plot(x_lens,y_plot,':k');
plot(x_lens,-y_plot,':k');

% Plot receiver plane
plot(xf*ones(1,2),lens,'-k','LineWidth',10)

% Plot boundaries
plot([xi,xf],lens(1)*ones(1,2),'--k')
plot([xi,xf],lens(2)*ones(1,2),'--k')

xlim([xi-0.5,xf+0.5]);
ylim([lens(1)-0.1,lens(2)+0.1]);

%% Plotting figure 2
figure(2); title('Lens'); hold on; movegui('east')

% Plot lens
plot(x_lens,y_plot,'-k');
plot(x_lens,-y_plot,'-k');

xlim([-0.5,0.5])
xlabel('$x$');
ylabel('$y$');

quiver(x_lens,y_plot,normal(:,1)',normal(:,2)')
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