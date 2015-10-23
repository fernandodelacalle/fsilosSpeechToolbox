clc; clear all; close all;

t_fin = 130;
t_ini = -10;
t_ini = -20;
time_step = 10;

frequencySize = 4*2;
frequencySize = 4*3;

[SE] = se3dNew(t_fin, t_ini, time_step, frequencySize);

SE = (SE - min(min(SE)))  ./ (max(max(SE)) -  min(min(SE))) ;

addZeros = ((t_fin + t_ini) / time_step) -1;
SE = [zeros(addZeros, frequencySize*2 + 1);  SE];

t = t_ini:time_step:t_fin;
bands = 3:-1/4:-3;
SE = SE';

figure(1);
surf(SE);
%xlabel('Time (mSeg)','fontsize',12);
ylabel('Frequency (Bands)','fontsize',12);




% figure(2);
% imagesc(SE);
% xlabel('Time (mSeg)');
% ylabel('Frequency (Bands)');


Nhood = ones(size(SE));
se_1 = strel(Nhood,SE);

