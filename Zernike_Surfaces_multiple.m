%% Zernike_Surfaces_multiple
% This code makes a plot of the surfaces of the zernike polynomials -
% together and single, and records a video to show a circular surface.

% Dependencies: CaptureFigVid.m
% Author:       Niels Buijssen 4561473
% Last updated: 28-04-2019

% Detailed description:
% Saves video to main folder. (not worth coding to save to other folder
% because it's only 1 video and don't expect to save more videos)

%% Settings
clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex');

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\Zernike Surfaces/';

Plot_multiple   = true;
Plot_single     = true;
Make_video      = false;

%% Input
n_radial = 1e2; % number of points plotted in radial directions
r = linspace(0,1,n_radial);
theta = linspace(0,2*pi,n_radial);
n_max = 4; % maximum n coef for zernike poly

%% Plotting

%% Plot multiple
if Plot_multiple == true
    
% Make grid for plotting
X = r.*sin(theta)';
Y = r.*cos(theta)';

% figure 1 topview

% angle for viewing plot
theta_view  = 0;  % corresponds to angle in x-y plane
phi_view    = 90; % corresponds to angle in x-z plane

figure('position', [0, 0, 1500, 1500]);
ind = 1;
for i = 1:(n_max+1)
    n = i - 1;
    m_arr = -n:2:n;
    for j = 1:length(m_arr)
        m = m_arr(j);
        subplot(n_max+1,n_max+1,ind)
            Z = Zer(n,m,r,theta);
            
            surf(X,Y,Z)
            colormap jet
            shading interp
            xlim([-1,1])
            ylim([-1,1])
            view(theta_view,phi_view)
            grid off
            axis off       
        ind = ind + 1;
    end
    ind = ind + n_max - length(m_arr)+1;
end

% Save figure
figure_name = 'Zernike_poly_topview';
filetype    = '.png';
print(figure(1), '-dpng', strcat(folder,figure_name,filetype))

% figure 2 top/sideview

% angle for viewing plot
theta_view  = 45;  % corresponds to angle in x-y plane
phi_view    = 45; % corresponds to angle in x-z plane

ind = 1;
for i = 1:(n_max+1)
    n = i - 1;
    m_arr = -n:2:n;
    for j = 1:length(m_arr)
        m = m_arr(j);
        subplot(n_max+1,n_max+1,ind)
            view(theta_view,phi_view)    
        ind = ind + 1;
    end
    ind = ind + n_max - length(m_arr)+1;
end
% Save figure
figure_name = 'Zernike_poly_top-sideview';
filetype    = '.png';
print(figure(1), '-dpng', strcat(folder,figure_name,filetype))

end

%% Plot single
if Plot_single == true
close all; figure(2)
    
% angle for viewing plot
theta_view  = 0;  % corresponds to angle in x-y plane
phi_view    = 90; % corresponds to angle in x-z plane


ind = 1;
for i = 1:(n_max+1)
    n = i - 1;
    m_arr = -n:2:n;
    for j = 1:length(m_arr)
        m = m_arr(j);
        Z = Zer(n,m,r,theta);
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
        ind = ind + 1;
        hold off
        
        % Save figure
        figure_name = strcat('Single_Zernike_',string(n),string(m));
        filetype    = '.png';
        print(figure(2), '-dpng', strcat(folder,figure_name,filetype))
    end
    ind = ind + n_max - length(m_arr)+1;
end

end

%% Make video
if Make_video == true
    
figure();
Z = -.3*Zer(2,0,r,theta)-0.3*Zer(3,-3,r,theta)+0.2*Zer(4,-2,r,theta);
X = r.*sin(theta)';
Y = r.*cos(theta)';
surf(X,Y,Z)
colormap jet
shading interp
xlim([-1,1])
ylim([-1,1])
zlim([0.8,1.3]);
xlabel('$x$');ylabel('$y$');zlabel('$z$');
view(45,45)
grid on
axis tight
OptionZ.FrameRate=40;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([0,45;90,45;180,45;270,45;362,90], 'Zernikeplot',OptionZ)

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