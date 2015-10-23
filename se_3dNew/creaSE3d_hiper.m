clc; clear all; close all;

%% Constants
t_fin = 200;
t_ini = -10;
%t_ini = -20;
time_step = 10;

frequencySize = 4*2;
frequencySize = 4*3;

SE = se3dAsymHyperboloid(t_fin, t_ini, time_step, frequencySize);
Nhood = ones(size(SE));
se_1 = strel(Nhood,SE);

figure(3)
surf(SE)
xlabel('Delta T')
ylabel('Delta F')
zlabel('Masking')