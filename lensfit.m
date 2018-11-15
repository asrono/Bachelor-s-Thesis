close all; clear all;
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
%% plot
figure(1)
% Plot optical axis
plot([0 l_OA],[0 0],'k--'); hold on

c = '.k';
for j = 1:l_N
    N_j = N_lens(:,:,j);
    X_j = X_lens(:,:,j);
    plot_point(N_j,c); 
    plot_point(X_j,c);
end

% Set graph options
title('SMS Method Lens Design')
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-3 3])
xlabel('Optical axis [m]')
ylabel('y [m]')
hold off

%% plot
figure(2);
plot3(N_3d(1,:),N_3d(2,:),N_3d(3,:),'.k'); hold on;
plot3(X_3d(1,:),X_3d(2,:),X_3d(3,:),'.k');
plot3(M_3d(1,:),M_3d(2,:),M_3d(3,:),'.k');
plot3(Y_3d(1,:),Y_3d(2,:),Y_3d(3,:),'.k');
xlim([-l_OA*0.02, l_OA*1.02])
ylim([-3 3]);
zlim([-3 3]);
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
grid on

%% functions
function plot_point(P,c)
    plot(P(1),P(2),c)
end

function in = extend_to_neg(in)
    in2 = in;
    in2(2,:,:) = -in(2,:,:);
    in(:,:,(end+1):(2*length(in))) = in2;
end