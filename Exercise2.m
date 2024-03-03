% Exercise 2
clear;close all;clc;
% 1. Define parameters
param.n = 100; % no. cells
param.D = 100; % depth
param.dz = param.D/param.n;

% 2. Set up grid
param.z =  param.dz/2:param.dz:param.D;

% 3. Set in. cond.
%P0 = normpdf(param.z,param.D/2,param.D/10)';
P0 = ones(100,1)*10;

% 4. Derivative function (time, phi)
param.u = 0; % current
param.d = 1; % Diffusion constant
%Ja = zeros(n+1,1);
%Jd = zeros(n+1,1);

% Set up odefun
% Parameters from paper
param.Iin = 350;
param.k = 15*10^(-12);
param.Kbg = 0.2;
param.l = 0.01;
%param.v = 0.04;
param.D = 5;
param.H = 30;
param.zm = 100;
param.zT = 20;
param.pmax = 0.04;
% 5. solve
tspan = 0:10000;
[t,P] = ode45(@(t,P)odefun2(t,P,param),tspan,P0);
P = P';

% 6. plot
figure(1)
image(P,'CDataMapping','scaled')
colorbar
xlabel('Time [Days]')
ylabel('Depth [meters]')
title('Conc. [10^3 cells/mL] over time')

%%
figure(2)
plot(P(:,50),1:100)
xlabel('Concentration')
ylabel('Depth [m]')
set ( gca, 'YDir', 'reverse' )