close all; clear all;
set(0,'defaulttextinterpreter','latex'); set(0,'defaultaxesfontsize',14);
%% Load data
savefile = 'lensini.mat';
load(savefile);
l_OA = ER(1,1,3); % length of optical axis
%% extend lens to negative y plane
N_lens = extend_to_neg(N_lens);
X_lens = extend_to_neg(X_lens);

% update length of N
l_N = length(N_lens);

% Input handling
N = zeros(2,l_N);
X = zeros(2,l_N);

N(:,:) = N_lens(:,1,:);
X(:,:) = X_lens(:,1,:);

% Functions
function in = extend_to_neg(in)
    in2 = in;
    in2(2,:,:) = -in(2,:,:);
    in(:,:,(end+1):(2*length(in))) = in2;
end