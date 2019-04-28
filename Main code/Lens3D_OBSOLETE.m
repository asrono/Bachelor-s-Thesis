%% Lens 3D [Obsolete]
% This code makes a plots to demonstrate taylor expansions.

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
savefile = 'lensini.mat';
load(savefile);
l_OA = ER(1,1,3); % length of optical axis
%% extend lens to negative y plane
N_lens = extend_to_neg(N_lens);
X_lens = extend_to_neg(X_lens);

% update length of N
l_N = length(N_lens);
%% 2d to 3d
N_3d = zeros(3,l_N);
X_3d = zeros(3,l_N);
M_3d = zeros(3,l_N);
Y_3d = zeros(3,l_N);

N_3d(1:2,:,:) = N_lens;
X_3d(1:2,:,:) = X_lens;

M_3d(1,:,:) = N_lens(1,:,:);
M_3d(3,:,:) = N_lens(2,:,:);
Y_3d(1,:,:) = X_lens(1,:,:);
Y_3d(3,:,:) = X_lens(2,:,:);

n = 16;
S1 = zeros(3,n*l_N);
for i = 1:n
theta = 2*pi*i/n;
A = [1, 0,           0;...
     0, cos(theta), -sin(theta);...
     0, sin(theta),  cos(theta)];
S1(:, (l_N*(i-1)+1):(l_N*i) ) = A*N_3d;
S2(:, (l_N*(i-1)+1):(l_N*i) ) = A*X_3d;
end

%% plot
figure(1); hold on;
plot3(S1(1,:),S1(2,:),S1(3,:),'.k','MarkerSize',5);
plot3(S2(1,:),S2(2,:),S2(3,:),'.k','MarkerSize',5);
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-3 3]);
zlim([-3 3]);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
set(gca,'ticklabelinterpreter','latex');
grid on
hold off;
%% functions
function plot_point(P,c)
    plot(P(1),P(2),c)
end

function in = extend_to_neg(in)
    in2 = in;
    in2(2,:,:) = -in(2,:,:);
    in(:,:,(end+1):(2*length(in))) = in2;
end