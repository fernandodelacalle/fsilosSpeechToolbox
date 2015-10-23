function [f, N_all] = frecuencyMaskingNew(N, centralBand)

NbansCoeff = 4; % 4 para 80 bandas 2 para 40

m = 60; 


a_pre = 30/NbansCoeff; 
a_pre = 4*NbansCoeff; %a_pre > N . With this interpretation is: 
                       %number of bands that takes to go to 0 from
                       %b_pre
b_pre = m - a_pre*centralBand;
N_pre = centralBand-N:centralBand;

f_pre = ellipse(a_pre, b_pre, N_pre); 
%a_pre*N_pre + b_pre;

a_pos = -8/NbansCoeff; 
a_pos = 3*NbansCoeff; %Same that with a_pre
b_pos = m - a_pos*centralBand;
N_pos = centralBand:centralBand+N;
f_pos = ellipse(a_pos,b_pre,N_pos); 
%a_pos*N_pos + b_pos;


f = [f_pre f_pos(2:end)];
N_all = [N_pre N_pos(2:end)];



end

