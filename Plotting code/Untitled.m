close all; clear all;

n_s = 1e3;
mu = 50;
sigma = 20;


x = normrnd(mu,sigma,[1,n_s]);

grades = 1 + 0.09*x;
grades2 = 1 + 0.1*x;

av1 = sum(grades)/n_s;
av2 = sum(grades2)/n_s;

figure(1);
hist(grades2,20,'g')

hold on;
hist(grades,20,'r');