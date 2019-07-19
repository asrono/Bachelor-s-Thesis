%% Results inverse problem
% Solves inverse problem and uses ray tracing and laplacian magic window to compare with the result.

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
%% Initialise systema

% At x = xi, we want to find a lens in zernike polynomials.
% At x = xf, we place a screen and observe intensity.

xi = 0; % arbitrary units
xf = 5; % arbitrary units

n_refractive = 1.6; % index of refraction

lens = [-1,1]; % position of bottom and top of lens
lens_middle = lens(1) + ( lens(2)-lens(1) );

rplot = linspace(0,1,1e3);
% intensity = 1 +0.2*cos(rplot*2*pi);
intensity = ones(1,1e3)+0.1*cos(rplot*2*pi);
% intensity = 0.8+0.5*normpdf(rplot,0,0.2)*sqrt(2*pi*0.2^2);

%% Find zernike coef
m = 0;
for n = 0:2:20
    a_vec_intensity(double2single_index(n,m)+1) = zernikecoef(n,m,rplot,intensity);
end
if false
clear a_vec_intensity intensity
loadfile = 'a_vec_Berry.mat';
load(strcat('C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Data\',loadfile),'a_vec_N');
a_vec = a_vec_N;    % Rename variable for this script
a_vec = a_vec/50;
intensity = 1 - (n_refractive - 1)*xf.*coef2surf(Zernike_laplacian_Janssen_vec(a_vec),rplot);
a_vec_intensity =  -(n_refractive - 1)*xf.*Zernike_laplacian_Janssen_vec(a_vec);
a_vec_intensity(1) = 1-0.01;
end


a_vec_lens = 1/(xf*(n_refractive-1))*Zernike_inverse_laplacian_Janssen_vec(a_vec_intensity);
a_vec_lens(5) = a_vec_lens(5) - 1/(8*(n_refractive-1)*xf);

lens_surface = coef2surf(a_vec_lens,rplot);
% for plotting use 2e2
number_of_rays = 1e5;

Der_ana = true; % true for analytical derivative calculation, false for numerical
plot_hist_3D = false; % plotting the 3D histogram (computationally heavy)

%% plotting
figure();
title('Intensity and fit'); hold on;
rplot2 = [fliplr(-rplot) rplot];
intensity2 = [fliplr(intensity) intensity];
intensity_fit = coef2surf(a_vec_intensity,rplot);
intensity_fit2 = [fliplr(intensity_fit) intensity_fit];
lens_surface2  = [fliplr(lens_surface) lens_surface];
plot(rplot2,intensity2,'-k');
plot(rplot2,intensity_fit2,'-r');
grid on;
xlabel('$r$');
ylabel('$I$');
legend({'Intensity','Zernike fit'},...
        'Interpreter','latex',...
        'Location', 'northwest');  



figure();
plot(rplot2,-lens_surface2,'-k')
title('Lens surface');
grid on;
xlabel('$r$');
ylabel('$h$');

I_lin = 1 - (n_refractive - 1)*xf.*coef2surf(Zernike_laplacian_Janssen_vec(-a_vec_lens),rplot);
I_lin = [fliplr(I_lin) I_lin];
I_nonlin = 1./ (1 + (n_refractive - 1)*xf.*coef2surf(Zernike_laplacian_Janssen_vec(-a_vec_lens),rplot));
I_nonlin = [fliplr(I_nonlin) I_nonlin];

figure(); hold on;
title('Intensity at $z = 5$')
plot(rplot2,I_lin,'-k');
plot(rplot2,I_nonlin,'--k');
grid on;
xlabel('$r$');
ylabel('$I$');
legend({'Linear','Non-linear'},...
        'Interpreter','latex',...
        'Location', 'northwest');  
    
%% Ray tracer
% for plotting
y_plot = linspace(lens(1),lens(2),number_of_rays);
x_lens = coef2surf(a_vec_lens,y_plot);
x_lens_deriv = coef2surf_deriv(a_vec_lens,y_plot);

%% Analytical differentiator
% calculate normal vectors at interface
dx = flipud(-coef2surf_deriv(a_vec_lens,y_plot)');
dy = ones(number_of_rays,1);

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

intensity = a.*2*xf + b;
% %% plotting figure 3
% figure(); title('Optical system'); hold on; movegui('center')
% 
% if number_of_rays < 2e2+1
% % Plot rays
% quiv2 = quiver(x_lens,y_plot,v_r(1,:),v_r(2,:),600);
% quiv2.Color = 'red';
% quiv2.LineWidth = 0.01;
% quiv2.MaxHeadSize = 0;
% 
% % Plot incoming rays
% quiv2 = quiver(x_lens,y_plot,-ones(1,number_of_rays),zeros(1,number_of_rays),60);
% quiv2.Color = 'red';
% quiv2.LineWidth = 0.01;
% quiv2.MaxHeadSize = 0;
% end
% % Plot optical axis
% plot([xi-2,xf],lens_middle*ones(1,2),'-k')
% 
% % % Plot lens % placeholder
% % plot(xi*ones(1,2),lens,':k')
% 
% % Plot receiver plane
% plot(xf*ones(1,2),lens,'-k','LineWidth',10)
% 
% % Plot boundaries
% plot([xi,xf],lens(1)*ones(1,2),'--k')
% plot([xi,xf],lens(2)*ones(1,2),'--k')
% 
% % Plot lens (real lens)
% plot(x_lens,y_plot,'-k');
% 
% xlabel('$x$')
% ylabel('$y$')
% xlim([xi-0.1,xf+0.5]);
% ylim([lens(1)-0.1,lens(2)+0.1]);
% 
% axis equal

figure();
title(strcat('Intensity at $z = ',string(xf),'$.')); hold on; movegui('west')

h1 = histogram(intensity,linspace(lens(1),lens(2),500),'EdgeColor','none');
xlim([lens(1),lens(2)])
ylabel('Counts / $I$');
xlabel('$r$');

% % Save figure
% figure_name = 'Ray_tracer1';
% filetype    = '.png';
% print(figure(1), '-dpng', strcat(folder,figure_name,filetype))
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
    c0 = 1./(r.^2.*(1-r.^2).^2);
    radial = c0 .* (...
        R(n  ,m  ,r).*(r.^4*(n+2)^2+2*r.^2*(m*(n+3)+n+2)+m^2)-...
        R(n+1,m+1,r).*(n+m+2).*r.*(r.^2*(2*n+6)+2*m+2)+...
        R(n+2,m+2,r).*(r.^2*(n+m+2)*(n+m+4)) );
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
