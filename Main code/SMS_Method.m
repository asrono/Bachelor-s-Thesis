%% SMS_Method
% The SMS (Simultaneous Multiple Surfaces) method can create points on the
% surfaces of a lens. This code gives a proof of concept and makes figures
% to show how the method works.

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 28-04-2019

% Detailed description:

%% Settings
close all; clear all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex'); 

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\SMS Method/';
%% Plotting options
plotHyp     = false;
plotLens    = true;
plotSym     = true;
plotRays    = false;
plotNormals = false;
SaveLens    = true;

%% Input parameters
% Etendue
U = 0.78; % [m^2]
% U < 2*[E_1,E_2] ( = 8 [m])

% Refraction of index lens
n = 1.6;

% angles for starting parabola
phi_E_start = 280.6;      % [degrees]
phi_E_end   = 360;      % [degrees]
phi_R_start = 0;        % [degrees]
phi_R_end   = 74;      % [degrees]
phi_1       = 281.3;    % [degrees]
phi_2       = 73.1    ; % [degrees]
n_phi       = 1e5;      % #points to plot

% Define emitter (E) and receiver (R) plane
E_1 = [0;-2];       % [m]
E_2 = [0; 2];       % [m]

l_OA = 30;          % [m]
R_1  = [l_OA; -1];  % [m]
R_2  = [l_OA;  1];  % [m]

%% Input handling
% Degrees to radians
phi_E_start = phi_E_start * pi/180;
phi_E_end   = phi_E_end   * pi/180;
phi_R_start = phi_R_start * pi/180;
phi_R_end   = phi_R_end   * pi/180;

phi_1   = phi_1   * pi/180;
phi_2   = phi_2   * pi/180;

phi_E = linspace(phi_E_start, phi_E_end, n_phi);
phi_R = linspace(phi_R_start, phi_R_end, n_phi);

% Emitter and receiver plane
ER = ones(1,2,4);
ER(:,:,1) = E_1;
ER(:,:,2) = E_2;
ER(:,:,3) = R_1;
ER(:,:,4) = R_2;

%% SMS Method
%% Step 1
h_E = hyp(U,E_1,E_2,phi_E);
h_R = hyp(U,R_1,R_2,phi_R);
N   = hyp(U,E_1,E_2,phi_1);
X   = hyp(U,R_1,R_2,phi_2);

%% Step 2
n_N = find_normal(E_1,N,X,n);
n_X = find_normal(R_1,X,N,n);

%% Step 3
% Find X_1 and N_1
X_1 = SMSs3(E_2,N,X,R_1,n,n_N);
N_1 = SMSs3(R_2,X,N,E_1,n,n_X);

n_X_1 = find_normal(R_1,X_1,N,n);
n_N_1 = find_normal(E_1,N_1,X,n);

%% Step 4
n_max = 100;
X_lens = zeros(2,1,n_max);
N_lens = zeros(2,1,n_max);
X_lens(:,:,1) = X;
N_lens(:,:,1) = N;
X_lens(:,:,2) = X_1;
N_lens(:,:,2) = N_1;

i = 2;

% Set n_max to prevent program for running too long. Second conditions
% stops the loop when the optical axis is reached.
while i < n_max && X_lens(2,1,i) >= 0
    X_i = X_lens(:,:,i);
    N_i = N_lens(:,:,i);
    [X_inc, N_inc] = SMSs4_layer(ER,N_i,X_i,n);
    X_lens(:,:,i+1) = X_inc;
    N_lens(:,:,i+1) = N_inc;
    i = i + 1;
end

% Remove extra points
N_lens = N_lens(:,:,1:i);
X_lens = X_lens(:,:,1:i);
%% Plotting
figure(1); % Step 1 and 2
% Plot transmitter and receiver planes
plot(E_1(1),E_1(2),'ro'); hold on
plot(E_2(1),E_2(2),'ro');
plot(R_1(1),R_1(2),'go');
plot(R_2(1),R_2(2),'go');

if plotHyp == true
% Plot hyperbolas
plot(h_E(1,:), h_E(2,:),'.k');
plot(h_R(1,:), h_R(2,:),'.b');
end

% Plot N and X points
if plotLens == true
c = '.k';
for j = 1:i
    N_j = N_lens(:,:,j);
    X_j = X_lens(:,:,j);
    plot_point(N_j,c);
    plot_point(X_j,c);
    if plotSym == true
        N_j(2) = -N_j(2);
        X_j(2) = -X_j(2);
        plot_point(N_j,c);
        plot_point(X_j,c);
    end
end
end

% Plot starting normals
if plotNormals == true
q_n_N = quiver(N(1),N(2),n_N(1),n_N(2));
q_n_X = quiver(X(1),X(2),n_X(1),n_X(2));

q_n_N.Color = 'black';
q_n_X.Color = 'blue';
q_n_N.LineWidth = 2;
q_n_X.LineWidth = 2;
q_n_N.MaxHeadSize = 0.8;
q_n_X.MaxHeadSize = 0.8;

q_n_N1 = quiver(N_1(1),N_1(2),n_N_1(1),n_N_1(2));
q_n_X1 = quiver(X_1(1),X_1(2),n_X_1(1),n_X_1(2));

q_n_N1.Color = 'black';
q_n_X1.Color = 'blue';
q_n_N1.LineWidth = 2;
q_n_X1.LineWidth = 2;
q_n_N1.MaxHeadSize = 0.8;
q_n_X1.MaxHeadSize = 0.8;
end

% Plot optical axis
plot([0 l_OA],[0 0],'k--')

% Plot rays
if plotRays == true
% Rays E_1
c = 'k';
plot_ray(E_1, N  , X  , R_1,c)
plot_ray(E_1, N_1, X  , R_2,c)

c = 'k';
plot_ray(E_2, N  , X_1, R_1,c)
end

% Set graph options
title('SMS Method Lens Design')
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-3 3])
xlabel('Optical axis [m]')
ylabel('y [m]')
hold off

% Save figure
figure_name = 'SMS_Overview';
filetype    = '.png';
print(figure(1), '-dpng', strcat(folder,figure_name,filetype))

%% Figure 2 starting points
figure(2);
E = [E_1 E_2];
R = [R_1 R_2];
% Plot transmitter and receiver planes
plot(E(1,:),E(2,:),'k','LineWidth',1); hold on
plot(R(1,:),R(2,:),'k','LineWidth',1);

% Plot hyperbolas
plot(h_E(1,:), h_E(2,:),'k','Linewidth',1);
plot(h_R(1,:), h_R(2,:),'k','Linewidth',1);

% Plot starting points
plot(N(1),N(2),'.k','Markersize',20)
plot(X(1),X(2),'.k','Markersize',20)

% Plot optical axis
plot([0 l_OA],[0 0],'k--','LineWidth',1)

% General settings
title('SMS method step 1')
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-1 1]*2.5)
axis off
hold off

% Annotations
text(E(1,1)-2.0,E(2,1),'$E_2$','FontSize',14);
text(E(1,2)-2.0,E(2,2),'$E_1$','FontSize',14);
text(R(1,1)+0.5,R(2,1),'$R_2$','FontSize',14);
text(R(1,2)+0.5,R(2,2),'$R_1$','FontSize',14);

text(5.9, 0.9, '$h_E$','FontSize',14);
text(l_OA-3.5, 0.9, '$h_R$','FontSize',14);

text(N(1)-2,N(2),'\textbf{N}','FontSize',14);
text(X(1)+0.5,X(2),'\textbf{X}','FontSize',14);

% Save figure
figure_name = 'SMS_step_1';
filetype    = '.png';
print(figure(2), '-dpng', strcat(folder,figure_name,filetype))

%% Figure 3 starting normals and rays
figure(3);
E = [E_1 E_2];
R = [R_1 R_2];
% Plot transmitter and receiver planes
plot(E(1,:),E(2,:),'k','LineWidth',1); hold on
plot(R(1,:),R(2,:),'k','LineWidth',1);

% Plot starting points
plot(N(1),N(2),'.k','Markersize',20)
plot(X(1),X(2),'.k','Markersize',20)

% Plot starting normals
pos = get(gca, 'Position');
annotation('arrow',[N(1),N(1)+ n_N(1)]/(l_OA*1.035),[1,1]/2.45+[N(2), N(2)+n_N(2)]/4.6);
annotation('arrow',[X(1),X(1)+ n_X(1)]/(l_OA*1.045),[1,1]/2.45+[X(2), X(2)+n_X(2)]/4.6);
text(N(1)-3,N(2)+0.5,'$\textbf{n}_N$','FontSize',14);
text(X(1)+1.5,X(2)+0.42,'$\textbf{n}_X$','FontSize',14);

% Plot ray 1
plot_ray_lw(E_1, N  , X  , R_1,'k',1)

% Plot optical axis
plot([0 l_OA],[0 0],'k--','LineWidth',1)

% General settings
title('SMS method step 2')
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-1 1]*2.5)
axis off
hold off

% Annotations
text(E(1,1)-2.0,E(2,1),'$E_2$','FontSize',14);
text(E(1,2)-2.0,E(2,2),'$E_1$','FontSize',14);
text(R(1,1)+0.5,R(2,1),'$R_2$','FontSize',14);
text(R(1,2)+0.5,R(2,2),'$R_1$','FontSize',14);

text(4, -0.8, '$r_1$','FontSize',14);
text(l_OA-2.5,-0.8, '$r_1$','FontSize',14);

text(N(1)-2,N(2),'\textbf{N}','FontSize',14);
text(X(1)+0.5,X(2),'\textbf{X}','FontSize',14);

% Save figure
figure_name = 'SMS_step_2';
filetype    = '.png';
print(figure(3), '-dpng', strcat(folder,figure_name,filetype))

%% Figure 4 second ray
figure(4);
E = [E_1 E_2];
R = [R_1 R_2];
% Plot transmitter and receiver planes
plot(E(1,:),E(2,:),'k','LineWidth',1); hold on
plot(R(1,:),R(2,:),'k','LineWidth',1);

% Plot starting points
plot(N(1),N(2),'.k','Markersize',20)
plot(X(1),X(2),'.k','Markersize',20)
plot(X_1(1),X_1(2),'.k','Markersize',20)

% Plot starting normals
pos = get(gca, 'Position');
annotation('arrow',[N(1),N(1)+ n_N(1)]/(l_OA*1.035),[1,1]/2.45+[N(2), N(2)+n_N(2)]/4.6);
annotation('arrow',[X(1),X(1)+ n_X(1)]/(l_OA*1.045),[1,1]/2.45+[X(2), X(2)+n_X(2)]/4.6);
text(N(1)-3,N(2)+0.5,'$\textbf{n}_N$','FontSize',14);
text(X(1)+1.5,X(2)+0.42,'$\textbf{n}_X$','FontSize',14);

% Plot ray 1
plot_ray_lw(E_1, N  , X  , R_1,'k',1)
% plot ray 2
plot_ray_lw(E_2, N, X_1 , R_1,'k',1)

% Plot optical axis
plot([0 l_OA],[0 0],'k--','LineWidth',1)

% General settings
title('SMS method step 3')
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-1 1]*2.5)
axis off
hold off

% Annotations
text(E(1,1)-2.0,E(2,1),'$E_2$','FontSize',14);
text(E(1,2)-2.0,E(2,2),'$E_1$','FontSize',14);
text(R(1,1)+0.5,R(2,1),'$R_2$','FontSize',14);
text(R(1,2)+0.5,R(2,2),'$R_1$','FontSize',14);

text(4, -0.8, '$r_1$','FontSize',14);
text(l_OA-2.5,-0.8, '$r_1$','FontSize',14);

text(4, 2.2, '$r_2$','FontSize',14);
text(l_OA-2.2,-0.2, '$r_2$','FontSize',14);

text(N(1)-2,N(2)+0.1,'\textbf{N}','FontSize',14);
text(X(1)+0.5,X(2),'\textbf{X}','FontSize',14)
text(X_1(1)-1,X_1(2)-0.3,'\textbf{X$_1$}','FontSize',14)

% Save figure
figure_name = 'SMS_step_3a';
filetype    = '.png';
print(figure(4), '-dpng', strcat(folder,figure_name,filetype))

%% Figure 5 second ray zoom
figure(5);
E = [E_1 E_2];
R = [R_1 R_2];
% Plot transmitter and receiver planes
plot(E(1,:),E(2,:),'k','LineWidth',1); hold on
plot(R(1,:),R(2,:),'k','LineWidth',1);

% Plot starting points
plot(N(1),N(2),'.k','Markersize',20)
plot(X(1),X(2),'.k','Markersize',20)
plot(X_1(1),X_1(2),'.k','Markersize',20)

% Plot starting normals
pos = get(gca, 'Position');
annotation('arrow',[N(1),N(1)+ n_N(1)]/(l_OA*1.55), [1,1]/2.45+[N(2), N(2)+n_N(2)]/5.5);
annotation('arrow',[X(1),X(1)+ n_X(1)]/(l_OA*1.105),[1,1]/2.45+[X(2), X(2)+n_X(2)]/5.9);
text(N(1)-3,N(2)+0.5,'$\textbf{n}_N$','FontSize',14);
text(X(1)+1.5,X(2)+0.42,'$\textbf{n}_X$','FontSize',14);

% Plot ray 1
plot_ray_lw(E_1, N  , X  , R_1,'k',1)
% plot ray 2
plot_ray_lw(E_2, N, X_1 , R_1,'k',1)

% Plot optical axis
plot([0 l_OA],[0 0],'k--','LineWidth',1)

% General settings
title('SMS method step 3')
xmargin = 1;
ymargin = 0.4;
xlim([N(1)-xmargin, X(1)+xmargin])
ylim([N(2)-ymargin, N(2)+ymargin/4])
axis off
hold off

% Annotations
text(19, 1.96, '$r_2$','FontSize',14);
text(20, 1.88, '$r_2$','FontSize',14);
text(21, 1.62, '$r_2$','FontSize',14);

text(19, 1.83, '$r_1$','FontSize',14);
text(20, 1.95, '$r_1$','FontSize',14);
text(21, 1.75, '$r_1$','FontSize',14);

text(4, 2.2, '$r_2$','FontSize',14);
text(l_OA-2.2,-0.2, '$r_2$','FontSize',14);

text(N(1)-0.2,N(2)+0.015,'\textbf{N}','FontSize',14);
text(X(1)+0.1,X(2),'\textbf{X}','FontSize',14)
text(X_1(1)+0.05,X_1(2),'\textbf{X$_1$}','FontSize',14)

annotation('textarrow',[0.25 0.15],[0.85 0.85],'interpreter','Latex','String','$E_1$','FontSize',14)
annotation('textarrow',[0.25 0.15],[0.24 0.12],'interpreter','Latex','String','$E_2$','FontSize',14)
annotation('textarrow',[0.8 0.9],[0.73 0.68],'interpreter','Latex','String','$R_1$','FontSize',14)
annotation('textarrow',[0.7 0.8],[0.25 0.09],'interpreter','Latex','String','$R_2$','FontSize',14)

% Save figure
figure_name = 'SMS_step_3b';
filetype    = '.png';
print(figure(5), '-dpng', strcat(folder,figure_name,filetype))

%% Figure 6 third ray zoom
figure(6);
E = [E_1 E_2];
R = [R_1 R_2];
% Plot transmitter and receiver planes
plot(E(1,:),E(2,:),'k','LineWidth',1); hold on
plot(R(1,:),R(2,:),'k','LineWidth',1);

% Plot starting points
plot(N(1),N(2),'.k','Markersize',20)
plot(X(1),X(2),'.k','Markersize',20)
plot(X_1(1),X_1(2),'.k','Markersize',20)
plot(N_1(1),N_1(2),'.k','Markersize',20)

% Plot ray 1
plot_ray_lw(E_1, N  , X  , R_1,':k',0.1)
% plot ray 2
plot_ray_lw(E_2, N, X_1 , R_1,':k',0.1)
% Plot ray 3
plot_ray_lw(E_1, N_1  , X  , R_2,'k',1)

% Plot optical axis
plot([0 l_OA],[0 0],'k--','LineWidth',1)

% General settings
title('SMS method step 3')
xmargin = 1;
ymargin = 0.4;
xlim([N(1)-xmargin, X(1)+xmargin])
ylim([N(2)-ymargin, N(2)+ymargin/4])
axis off
hold off

% Annotations
text(19, 1.77, '$r_3$','FontSize',14);
text(20, 1.92, '$r_3$','FontSize',14);
text(21, 1.88, '$r_3$','FontSize',14);

text(N(1)-0.2,N(2)+0.015,'\textbf{N}','FontSize',14);
text(N_1(1)-0.2,N_1(2)+0.015,'\textbf{N$_1$}','FontSize',14);
text(X(1)+0.1,X(2),'\textbf{X}','FontSize',14)
text(X_1(1)+0.05,X_1(2),'\textbf{X$_1$}','FontSize',14)

annotation('textarrow',[0.25 0.15],[0.85 0.85],'interpreter','Latex','String','$E_1$','FontSize',14)
annotation('textarrow',[0.25 0.15],[0.24 0.12],'interpreter','Latex','String','$E_2$','FontSize',14)
annotation('textarrow',[0.8 0.9],[0.73 0.68],'interpreter','Latex','String','$R_1$','FontSize',14)
annotation('textarrow',[0.7 0.8],[0.25 0.09],'interpreter','Latex','String','$R_2$','FontSize',14)

% Save figure
figure_name = 'SMS_step_3c';
filetype    = '.png';
print(figure(6), '-dpng', strcat(folder,figure_name,filetype))

%% Figure 7 fourth and fifth ray zoom
figure(7); close; figure(7);
E = [E_1 E_2];
R = [R_1 R_2];
% Plot transmitter and receiver planes
plot(E(1,:),E(2,:),'k','LineWidth',1); hold on
plot(R(1,:),R(2,:),'k','LineWidth',1);

% Plot starting points
plot(N(1),N(2),'.k','Markersize',20)
plot(X(1),X(2),'.k','Markersize',20)
plot(X_1(1),X_1(2),'.k','Markersize',20)
plot(N_1(1),N_1(2),'.k','Markersize',20)
plot(X_lens(1,1,3),X_lens(2,1,3),'.k','Markersize',20);
plot(N_lens(1,1,3),N_lens(2,1,3),'.k','Markersize',20);

% Plot ray 1
plot_ray_lw(E_1, N  , X  , R_1,':k',0.1)
% plot ray 2
plot_ray_lw(E_2, N, X_1 , R_1,':k',0.1)
% Plot ray 3
plot_ray_lw(E_1, N_1  , X  , R_2,':k',0.1)
% plot ray 4
plot_ray_lw(E_2, N_1, X_lens(:,:,3) , R_1,'k',1)
% plot ray 5
plot_ray_lw(E_1, N_lens(:,:,3), X_lens(:,:,2) , R_2,'k',1)


% Plot optical axis
plot([0 l_OA],[0 0],'k--','LineWidth',1)

% General settings
title('SMS method step 4')
xmargin = 1;
ymargin = 0.4;
xlim([N(1)-xmargin, X(1)+xmargin])
ylim([N(2)-ymargin, N(2)+ymargin/4])
axis off
hold off

% Annotations
text(19, 1.89, '$r_4$','FontSize',14);
text(19.86, 1.855, '$r_4$','FontSize',14);
text(21, 1.61, '$r_4$','FontSize',14);

text(19, 1.69, '$r_5$','FontSize',14);
text(19.86, 1.78, '$r_5$','FontSize',14);
text(21, 1.81, '$r_5$','FontSize',14);
text(N(1)-0.2,N(2)+0.015,'\textbf{N}','FontSize',14);
text(N_1(1)-0.2,N_1(2)+0.015,'\textbf{N$_1$}','FontSize',14);
text(N_lens(1,1,3)-0.2,N_lens(2,1,3)+0.015,'\textbf{N$_2$}','FontSize',14);
text(X(1)+0.1,X(2),'\textbf{X}','FontSize',14)
text(X_1(1)+0.05,X_1(2)+0.01,'\textbf{X$_1$}','FontSize',14)
text(X_lens(1,1,3)+0.05,X_lens(2,1,3)+0.01,'\textbf{X$_2$}','FontSize',14)

annotation('textarrow',[0.25 0.15],[0.85 0.85],'interpreter','Latex','String','$E_1$','FontSize',14)
annotation('textarrow',[0.25 0.15],[0.24 0.12],'interpreter','Latex','String','$E_2$','FontSize',14)
annotation('textarrow',[0.8 0.9],[0.73 0.68],'interpreter','Latex','String','$R_1$','FontSize',14)
annotation('textarrow',[0.7 0.8],[0.25 0.09],'interpreter','Latex','String','$R_2$','FontSize',14)

% Save figure
figure_name = 'SMS_step_4a';
filetype    = '.png';
print(figure(7), '-dpng', strcat(folder,figure_name,filetype))

%% Figure 8 whole lens
figure(8);
E = [E_1 E_2];
R = [R_1 R_2];

% Plot transmitter and receiver planes
plot(E(1,:),E(2,:),'k','LineWidth',1); hold on
plot(R(1,:),R(2,:),'k','LineWidth',1);


% Plot points
plot(N(1),N(2),'.k','Markersize',20)
plot(X(1),X(2),'.k','Markersize',20)

for j = 1:(i-1)
    N_j = N_lens(:,:,j);
    X_j = X_lens(:,:,j);
    plot(N_j(1),N_j(2),'.k','Markersize',20)
    plot(X_j(1),X_j(2),'.k','Markersize',20)
    plot(N_j(1),-N_j(2),'ok','Markersize',5)
    plot(X_j(1),-X_j(2),'ok','Markersize',5)
end

% Plot optical axis
plot([0 l_OA],[0 0],'k--','LineWidth',1)

% General settings
title('SMS method complete')
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-1 1]*2.5)
axis off
hold off

% Annotations
text(E(1,1)-2.0,E(2,1),'$E_2$','FontSize',14);
text(E(1,2)-2.0,E(2,2),'$E_1$','FontSize',14);
text(R(1,1)+0.5,R(2,1),'$R_2$','FontSize',14);
text(R(1,2)+0.5,R(2,2),'$R_1$','FontSize',14);

text(N(1)-2,N(2),'\textbf{N}','FontSize',14);
text(X(1)+0.5,X(2),'\textbf{X}','FontSize',14);
text(N(1)-2,-N(2),'\textbf{M}','FontSize',14);
text(X(1)+0.5,-X(2),'\textbf{Y}','FontSize',14);

% Save figure
figure_name = 'SMS_step_4b';
filetype    = '.png';
print(figure(8), '-dpng', strcat(folder,figure_name,filetype))

%% Save variables
if SaveLens
    savefile = 'lensini.mat';
    save(savefile,'N_lens', 'X_lens', 'ER');
end
%% Functions
function plot_ray(E,G,F,R,c)
    plot([E(1) G(1) F(1) R(1)],[E(2) G(2) F(2) R(2)],c)
end

function plot_ray_lw(E,G,F,R,c,lw)
    plot([E(1) G(1) F(1) R(1)],[E(2) G(2) F(2) R(2)],c,'LineWidth',lw)
end

function plot_ray_lw_clr(E,G,F,R,c,lw)
    plot([E(1) G(1) F(1) R(1)],[E(2) G(2) F(2) R(2)],'Color',c,'LineWidth',lw)
end

function plot_point(P,c)
    plot(P(1),P(2),c)
end
function h = hyp(U,P_1,P_2,phi)
    alpha = angh(P_2-P_1);
    h = P_1 + ( (U/2)^2  - norm(P_1-P_2)^2 ) ./ (U - 2 * norm(P_1-P_2) * cos(phi) ) .* [cos(phi+alpha); sin(phi+alpha)];
end

function a = ang(v,u)
    a = acos(dot(v,u)/(norm(v)*norm(u)));
end

function a = angp(v,u)
    if numel(v) ~= 2 || numel(u) ~= 2
        error('Angle can only be determined in 2D')
    end
    
    if u(1)*v(2)-u(2)*v(1) >= 0
        a = ang(v,u);
    else
        a = 2*pi - ang(v,u);
    end
end

function a = angh(v)
    a = angp(v,[1 0]);
end

function normal = find_normal(E,N,X,n)
    A = N - E;
    B = X - N;
    C = A/norm(A)-n*B/norm(B);
    normal = C/norm(C);
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
    v_r = v_r/norm(v_r);
end

function Q = find_new_point(G,F,P,v_r,n)
    S = n*norm(G-F)+norm(F-P);
    C_1 = n*S + dot(G-P,v_r);
    C_2 = S^2-norm(G-P)^2;
    Q   = G + v_r*(C_1 - sqrt(C_2*(1-n^2)+C_1^2))/(n^2-1);
end

function F_i = SMSs3(E,G,F,R,n,n_G)
    v_i = (G-E)/norm(G-E);
    v_r = find_reflected_ray(v_i,1,n,n_G);
	F_i = find_new_point(G,F,R,v_r,n);  
end

function X_inc = SMSs4(ER,N_i,X_i,n,right)
    if right == 1
        E_1 = ER(:,:,1)';
        E_2 = ER(:,:,2)';
        R_1 = ER(:,:,3)';
    else
        E_1 = ER(:,:,3)';
        E_2 = ER(:,:,4)';
        R_1 = ER(:,:,1)';
        [N_i, X_i] = deal(X_i,N_i);
    end
    
    n_N = find_normal(E_1,N_i,X_i,n);
    X_inc = SMSs3(E_2,N_i,X_i,R_1,n,n_N);
end

function [X_inc, N_inc] = SMSs4_layer(ER,N_i,X_i,n)
    X_inc = SMSs4(ER,N_i,X_i,n,1);
    N_inc = SMSs4(ER,N_i,X_i,n,0);
end