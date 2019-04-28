close all; clear all;
set(0,'defaulttextinterpreter','latex'); set(0,'defaultaxesfontsize',14);
%% Load data
savefile = 'lensini.mat';
load(savefile);
l_OA = ER(1,1,3); % length of optical axis
% %% extend lens to negative y plane
% N_lens = extend_to_neg(N_lens);
% X_lens = extend_to_neg(X_lens);
% 
% update length of N
l_N = length(N_lens);

% figure(1);
% plot_points(N_lens,'ok')

f = zeros(1,l_N);
r = f;

f(1,:) = N_lens(1,:,:);
r_0(1,:) = N_lens(2,:,:);
r(1,:) = N_lens(2,:,:)/max(r_0);

r = r(1:end-1);
f = f(1:end-1);
l_N = l_N-1;
figure(2);
plot(r,f,'ko')
hold on
r_int = linspace(0,1,1e3);
f_int = interp1(r,f,r_int,'spline');

r = r_int;
f = f_int;

m = 0;
for n = 0:2:6
    double2single_index(n,m);
    a_vec(double2single_index(n,m)+1) = zernikecoef(n,m,r,f);
end

r_prime = 0:0.01:max(r);
for j = 0:(length(a_vec)-1)
    [n,m] = single2double_index(j);
    rad_mat(j+1,:) = a_vec(j+1)*R(n,m,r_prime);
end
rad = sum(rad_mat,1);

plot(r_prime,rad,'r')
%%
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