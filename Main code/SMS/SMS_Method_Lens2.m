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
plotHyp     = true;
plotLens    = true;
plotSym     = true;
plotRays    = false;
plotNormals = false;
SaveLens    = true;

%% Input parameters
% Etendue
U = 0.3; % [m^2]
% U < 2*[E_1,E_2] ( = 8 [m])

% Refraction of index lens
n = 1.6;

% angles for starting parabola
phi_E_start = 322;        % [degrees]
phi_E_end   = 321;      % [degrees]
phi_R_start = 83;       % [degrees]
phi_R_end   = 85;      % [degrees]
phi_1       = 321.1;  % [degrees]
phi_2       = 83.89;   % [degrees]
n_phi       = 1e5;      % #points to plot

% Define emitter (E) and receiver (R) plane
E_1 = [0;-0.1];       % [m]
E_2 = [0; 0.1];       % [m]

l_OA = 30;             % [m]
R_1  = [l_OA; -2];  % [m]
R_2  = [l_OA;  2];  % [m]

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
n_max = 1000;
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

%% Save variables
if SaveLens
    savefile = 'lensExample2.mat';
    save(strcat('C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Data\',savefile),'N_lens', 'X_lens', 'ER');
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
