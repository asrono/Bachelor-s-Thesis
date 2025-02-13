%% Results Direct Problem
% Solves the direct problem via Laplacian magic window and ray tracing.

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 21-06-2019

% Detailed description:

%% Settings
close all; clear all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex'); 

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\Ray Tracer Zernike\';

%% General Settings

% choose 1
input_lens = false;
create_lens = true;

%% Check validity settings
if ~xor(input_lens == true, create_lens == true)
    error('Set either Input_lens or Create_lens to true')
end

n_fit = 1e3;
r = linspace(0,1,n_fit);
%% Input lens
if input_lens == true
% load coefficients for a lens
loadfile = 'a_vec_Berry.mat';
load(strcat('C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Data\',loadfile),'a_vec_N');
a_vec = a_vec_N;    % Rename variable for this script
a_vec(1) = 0;       % Set height of lens to optical start
a_vec = a_vec/50;
clear a_vec_N
end

%% Create lens
if create_lens == true
% h = 1/200*normpdf(r,0,0.3);
h = 1/2000*cos(4*pi*r);

% Find zernike coef
m = 0;
n_max = 20;
for n = 0:2:n_max
    a_vec(double2single_index(n,m)+1) = zernikecoef(n,m,r,h);
end

% plot lens and fit
figure(); hold on;
plot(r,h,'-k','linewidth',2)
plot(r,coef2surf(a_vec,r),'-r','linewidth',1.5)
end

%% Set up for intensity
zf = 30; % Final z value
z = linspace(0,zf,1e3);
n_ref = 1.6; % index of refraction

% calculate intensity for different z
% Berry linear
I_lin = 1 - (n_ref - 1)*z'.*coef2surf(Zernike_laplacian_Janssen_vec(a_vec),r);
I_lin = [fliplr(I_lin) I_lin];

% Berry non linear
I_nonlin = 1./ (1 + (n_ref - 1)*z'.*coef2surf(Zernike_laplacian_Janssen_vec(a_vec),r));
I_nonlin = [fliplr(I_nonlin) I_nonlin];

% update r
r2 = [fliplr(-r) r];
%% Plot intensity
figure();
subplot(2,2,1)
surf(r2,z,I_lin);
shading interp
view(-35,25);
xlabel('$r$');
ylabel('$z$');
zlabel('$I$');
title('Linear approximation')
zlim([0.5,2]);

subplot(2,2,2)
plot(r2,I_lin(end,:),'-k'); hold on;
plot(r2,I_nonlin(end,:),'--k')
xlabel('$r$');
ylabel('$I$');
grid on;
title('Intensity at $z = 5$')
legend({'Linear','Non-linear'},...
        'Interpreter','latex',...
        'Location', 'northwest');  
subplot(2,2,3)
surf(r2,z,I_nonlin);
shading interp
view(-35,25);
xlabel('$r$');
ylabel('$z$');
zlabel('$I$');
title('Non-linear approximation')
zlim([0.5,2])

subplot(2,2,4)
plot(r2,abs(I_nonlin(end,:)-I_lin(end,:)),'-k');
xlabel('$r$');
ylabel('$I$');
grid on;
title('Diff. intensity at $z = 5$')

% figure();
% plot(z,sum(I_lin,2)/2e3)
%% Initialise system

% At x = xi, we place a lens with known zernike polynomials.
% At x = xf, we place a screen and observe intensity.

xi = 0; % arbitrary units
xf = zf; % arbitrary units

lens = [-1,1]; % position of bottom and top of lens
lens_middle = lens(1) + ( lens(2)-lens(1) );


%% Settings
% for plotting use 2e2
number_of_rays = 1e5;

Der_ana = true ; % true for analytical derivative calculation, false for numerical
plot_hist_3D = true; % plotting the 3D histogram (computationally heavy)

%% Ray tracer
% for plotting
y_plot = linspace(lens(1),lens(2),number_of_rays);
x_lens = coef2surf(a_vec,y_plot);
x_lens_deriv = coef2surf_deriv(a_vec,y_plot);

rays_start = linspace(lens(2),lens(1),number_of_rays)';

% %% Numerical differentiator
% if Der_ana == false
% % calculate normal vectors at interface
% dx = zeros(number_of_rays,1);
% dx(1) =   x_lens(2)     - x_lens(1);
% dx(end) = x_lens(end) - x_lens(end-1);
% for i = 2: (number_of_rays-1)
%    dx(i) = x_lens(i+1)-x_lens(i-1);
% end
% dy = ones(number_of_rays,1) * 2 *( lens(1) - lens(2) ) / (number_of_rays - 1);
% dy(1) = dy(1)/2;
% dy(end) = dy(end)/2;
% dy = -dy;
% end
%% Analytical differentiator
if Der_ana == true
% calculate normal vectors at interface
dx = flipud(-coef2surf_deriv(-a_vec,y_plot)');
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
%%
%% plotting figure 1
figure(); title('Optical system'); hold on; movegui('center')

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
xlim([xi-0.1,xf]);
ylim([lens(1)-0.1,lens(2)+0.1]);
% axis equal
% 
% % Save figure
% figure_name = 'Ray_tracer1';
% filetype    = '.png';
% print(figure(1), '-dpng', strcat(folder,figure_name,filetype))
% 
% %% Plotting figure 2
% figure(2); title('Lens'); hold on; movegui('east')
% 
% % Plot lens
% plot(x_lens,y_plot,'-k');
% plot(x_lens_deriv,y_plot,'-r');
% 
% if number_of_rays < 2e2+1
% % Plot normals on lens
% quiv = quiver(x_lens,y_plot,normal(1,:),normal(2,:));
% quiv.Color = 'black';
% quiv.LineWidth = 1;
% quiv.MaxHeadSize = 0.5;
% end
% 
% xlim([-0.5,0.5])
% xlabel('$x$');
% ylabel('$y$');
% 
% grid minor
% axis equal
% 
%% Plotting figure 3
figure(); title(strcat('Intensity at $z = ',string(zf),'$.')); hold on; movegui('west')

h1 = histogram(intensity,linspace(lens(1),lens(2),500),'EdgeColor','none');
xlim([lens(1),lens(2)])
ylabel('Counts / $I$')
xlabel('$r$')
% Save figure
figure_name = 'ray_tracer3';
filetype    = '.png';
print(figure(3), '-dpng', strcat(folder,figure_name,filetype))

% figure();
% midpoints = conv(h1.BinEdges, [0.5 0.5], 'valid');
% counts = h1.BinCounts;
% plot(midpoints,counts,'or'); hold on 
% % find coefficients
% m = 0;
% midpoints2 = midpoints(round(end/2):end);
% counts2 = counts(round(end/2):end)+ 0.00001*rand(1,length(counts(round(end/2):end)));
% plot(midpoints2,counts2,'ok')
% % set amount of interpolated points
% n_interpolation = 1e5; % recommended > 1e3
% % get interpolated points
% [midpoints3, counts3] = ...
%     interpol(counts2, midpoints2, 0,1,n_interpolation);
% plot(midpoints3,counts3,'.g')
% for n = 0:2:30
%     a_vec_2(double2single_index(n,m)+1) = zernikecoef(n,m,midpoints2,counts2); %#ok<SAGROW>
% end
% plot(midpoints3,coef2surf(a_vec_2,midpoints3),'-r')
% figure(3);


%% Plotting figure 4
if plot_hist_3D == true
figure();
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
xlabel('$r$')
ylabel('$z$')
zlabel('Counts')
end

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
