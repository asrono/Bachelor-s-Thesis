%% Zernike_Surfaces_Single
% This code makes plots of the surfaces of single zernike polynomials.

% Dependencies: None
% Author:       Niels Buijssen 4561473
% Last updated: 28-04-2019

% Detailed description:

%% Settings
clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaultAxesTickLabelInterpreter','latex');

folder = 'C:\Users\Buijssen\Documents\GitHub\Bachelor-s-Thesis\Figures\Zernike Surfaces/';

%% Input
n_radial = 1e2; % number of points plotted in radial directions
r = linspace(0,1,n_radial);
theta = linspace(0,2*pi,n_radial);

%% Input handling
X = r.*sin(theta)';
Y = r.*cos(theta)';

%% Functions
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