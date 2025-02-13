%% Ray tracer Zernike
% Creates ray tracer for Zernike purposes and laplacian magic window

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 14-05-2019

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
% for plotting use 2e2
number_of_rays = 1e4;

Der_ana = true; % true for analytical derivative calculation, false for numerical
plot_hist_3D = false; % plotting the 3D histogram (computationally heavy)

%% Input
% load coefficients for a lens
loadfile = 'a_vec_Berry.mat';
load(strcat('C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Data\',loadfile),'a_vec_N');
a_vec = a_vec_N;    % Rename variable for this script
a_vec(1) = xi;       % Set height of lens to optical start
a_vec = a_vec/2;
% a_vec = zeros(1,25);
% a_vec(12) = 0.01;
% a_vec(13) = 0.02;
% % a_vec(25) = 0.01;
clear a_vec_N

%% Ray tracer
% for plotting
y_plot = linspace(lens(1),lens(2),number_of_rays);
x_lens = coef2surf(a_vec,y_plot);
x_lens_deriv = coef2surf_deriv(a_vec,y_plot);

rays_start = linspace(lens(2),lens(1),number_of_rays)';

%% Numerical differentiator
if Der_ana == false
% calculate normal vectors at interface
dx = zeros(number_of_rays,1);
dx(1) =   x_lens(2)     - x_lens(1);
dx(end) = x_lens(end) - x_lens(end-1);
for i = 2: (number_of_rays-1)
   dx(i) = x_lens(i+1)-x_lens(i-1);
end
dy = ones(number_of_rays,1) * 2 *( lens(1) - lens(2) ) / (number_of_rays - 1);
dy(1) = dy(1)/2;
dy(end) = dy(end)/2;
dy = -dy;
end
%% Analytical differentiator
if Der_ana == true
% calculate normal vectors at interface
dx = flipud(-coef2surf_deriv(a_vec,y_plot)');
dy = ones(number_of_rays,1);
end

normal = [dy,-dx];
norm = sqrt( normal(:,1).^2 + normal(:,2).^2 );
normal = (normal ./ norm)';

% calculate rays after transmission
n_i = 1.6;
n_r = 1.0;
v_i = [1,0];
v_r = zeros(2,number_of_rays);
for i = 1:number_of_rays
    v_r(:,i) = find_reflected_ray(v_i,n_i,n_r,normal(:,i)');
end


a = v_r(2,:)./v_r(1,:);
b = y_plot - a.*x_lens;

intensity = a.*xf + b;
%% plotting figure 1
figure(1); title('Optical system'); hold on; movegui('center')

if number_of_rays < 2e2+1
% Plot rays
quiv2 = quiver(x_lens,y_plot,v_r(1,:),v_r(2,:),600);
quiv2.Color = 'red';
quiv2.LineWidth = 0.01;
quiv2.MaxHeadSize = 0;

% Plot incoming rays
quiv2 = quiver(x_lens,y_plot,-ones(1,number_of_rays),zeros(1,number_of_rays),60);
quiv2.Color = 'red';
quiv2.LineWidth = 0.01;
quiv2.MaxHeadSize = 0;
end
% Plot optical axis
plot([xi-2,xf],lens_middle*ones(1,2),'-k')

% % Plot lens % placeholder
% plot(xi*ones(1,2),lens,':k')

% Plot receiver plane
plot(xf*ones(1,2),lens,'-k','LineWidth',10)

% Plot boundaries
plot([xi,xf],lens(1)*ones(1,2),'--k')
plot([xi,xf],lens(2)*ones(1,2),'--k')

% Plot lens (real lens)
plot(x_lens,y_plot,'-k');

xlabel('$x$')
ylabel('$y$')
xlim([xi-0.1,xf+0.5]);
ylim([lens(1)-0.1,lens(2)+0.1]);
% axis equal

% Save figure
figure_name = 'Ray_tracer1';
filetype    = '.png';
print(figure(1), '-dpng', strcat(folder,figure_name,filetype))

%% Plotting figure 2
figure(2); title('Lens'); hold on; movegui('east')

% Plot lens
plot(x_lens,y_plot,'-k');
plot(x_lens_deriv,y_plot,'-r');

if number_of_rays < 2e2+1
% Plot normals on lens
quiv = quiver(x_lens,y_plot,normal(1,:),normal(2,:));
quiv.Color = 'black';
quiv.LineWidth = 1;
quiv.MaxHeadSize = 0.5;
end

xlim([-0.5,0.5])
xlabel('$x$');
ylabel('$y$');

grid minor
axis equal

%% Plotting figure 3
figure(3); title(strcat('Intensity at $x = x_f = ',string(xf),'$.')); hold on; movegui('west')

h1 = histogram(intensity,linspace(lens(1),lens(2),100));
xlim([lens(1),lens(2)])
ylabel('Counts (\#rays)')
xlabel('$x$')
% Save figure
figure_name = 'ray_tracer3';
filetype    = '.png';
print(figure(3), '-dpng', strcat(folder,figure_name,filetype))

figure(5);
midpoints = conv(h1.BinEdges, [0.5 0.5], 'valid');
counts = h1.BinCounts;
plot(midpoints,counts,'or'); hold on 
% find coefficients
m = 0;
midpoints2 = midpoints(round(end/2):end);
counts2 = counts(round(end/2):end)+ 0.00001*rand(1,length(counts(round(end/2):end)));
plot(midpoints2,counts2,'ok')
% set amount of interpolated points
n_interpolation = 1e5; % recommended > 1e3
% get interpolated points
[midpoints3, counts3] = ...
    interpol(counts2, midpoints2, 0,1,n_interpolation);
plot(midpoints3,counts3,'.g')
for n = 0:2:20
    a_vec_2(double2single_index(n,m)+1) = zernikecoef(n,m,midpoints2,counts2); %#ok<SAGROW>
end
plot(midpoints3,coef2surf(a_vec_2,midpoints3),'-r')
figure(3);


%% Plotting figure 4
if plot_hist_3D == true
figure(4);
number_x = 5e2;
x = linspace(0,xf,number_x)';
t = a.*x + b;

binNo = 3e2;                     % number of bins for the histogram
bins  = zeros(size(t,1),binNo); % preallocation 
% 'FOR' loop to get the histogram values for each photo
for i=1:size(t,1)
h = histogram(t(i,:),linspace(lens(1),lens(2),binNo+1));
bins(i,:) = h.Values;
end

b = bar3(bins);

for i = 1:size(bins,2)
    cdata = get(b(i),'cdata');
    k = 1;
    for j = 0:6:(6*size(bins,1)-6)
        cdata(j+1:j+6,:) = bins(k,i);
        k = k+1;
    end
    set(b(i),'cdata',cdata);
end
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
colormap jet
shading interp
view(180,90)
axis tight
xlabel('$x$')
ylabel('$y$')
zlabel('Counts')
end

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

function radial = R_deriv2(n,m,r)
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
