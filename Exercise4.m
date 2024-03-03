% Exercise 4
clear;close all;clc;
% 1. Define parameters
param.n = 100; % no. cells
param.D = 100; % depth
param.dz = param.D/param.n; % grid size
param.z =  param.dz/2:param.dz:param.D; % grid

param.u = 0.04; % plankton current (gravitational flux) / sinking velocity (m/h)
param.d = 1; % Diffusion constant (cm2/s)
param.Iin = 400; % Light surface input (micromol photons/m2 s)
param.k = 0.01*6*10^(-10); % light absorbtion coefficient of plankton (m2/cell)
param.Kbg = 0.045; % background turbidity (1/m)
param.HI = 20; % Half saturation constant of light limited growth (micromol photons/m2 s)
param.HN = 0.015; % Half saturation constant of nutrient limited growth (mmol nutrient/m3)
param.pmax = 0.04; % maximum growth rate (both plankton and nutrients) (1/h)
param.y = 1*10^(-9); % Nutrient content of phytoplankton (mmol nutrient/cell)
param.l = 0.01; % specific loss rate (plankton mortality) (1/h)
param.ND = 5; % Nutrient conc. at bottom (mmol nutrient/m3)
param.gamma = 0.5/24; % Grazing (m3 (mmol N)^(-1) 1/h)
param.r = 0.15/24; % Detritus remineralization rate (1/h)
param.w = 10/24; % detritus sinking velocity (1/h)

% 5. solve
P0 = ones(param.n,1)*10; % Initial conditions:
N0 = ones(param.n,1)*5;
D0 = ones(param.n,1);
Y0 = [P0;N0;D0];
tspan = 0:1:5000;
[t,Y] = ode45(@(t,Y)odefun4(t,Y,param),tspan,Y0);
[~, pI,pN] = odefun4(t(end),Y(end,:),param); % Finding the limiting factors (light or minerals)

Y = Y';
P = Y(1:param.n,:);
N = Y(param.n+1:(2*param.n),:);
D = Y((2*param.n+1):end,:);

% 6. plot
% Plankton
figure(1)
image(P,'CDataMapping','scaled')
colorbar
xlabel('Time [100 hours]')
ylabel('Depth [meters]')
title('Phytoplankton conc. [cells/m3] over time')

% Nutrients
figure(2)
image(N,'CDataMapping','scaled')
colorbar
xlabel('Time [hours]')
ylabel('Depth [meters]')
title('Nutrient conc. [mmol nutrient/m3]')

% Detritus
figure(3)
image(D,'CDataMapping','scaled')
colorbar
xlabel('Time [hours]')
ylabel('Depth [meters]')
title('Detritus conc. [mmol detritus/m3]')

% Limiting factors
figure(4)
plot(pI,param.z)
set(gca, 'YDir','reverse')
hold on
plot(pN,param.z)
set(gca, 'YDir','reverse')
hold off
legend('Light','Minerals')
title('Limiting factors')
ylabel('Depth [meters]')
