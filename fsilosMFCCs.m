%[c] = fsilosMFCCs(Signal, samplingRate, SS_flag, MF_flag, se)
function [c] = fsilosMFCCs(Signal, samplingRate, SS_flag, MF_flag, se)

if nargin < 2;   error('s and fs requiered'); return;  end   
if nargin < 3;   SS_flag = 0;   end
if nargin < 4;   MF_flag = 0; se = 0;  end 
if nargin < 5;   
    if MF_flag
        load('SEs/elementoFran.mat');
        SE = SE';
        Nhood = ones(size(SE));
        se = strel(Nhood,SE);
    else
        se = 0;
    end
end

% Tipo de filtros
filterbank = 'mfcc';
% Numero de bandas
nbands = 40;  
[fb,aspc] = fsilosMelfCC(Signal, samplingRate, filterbank, MF_flag,  se,  SS_flag, 'numcep', 13, 'maxfreq', samplingRate/2, 'nbands', nbands, 'fbtype', 'fcmel', 'dcttype', 2, 'preemph', 0.97, 'dither', 1);


% Los dos pasos a continuación: El primero se hace al final, y el segundo
% es innecesario, pues se normaliza en el script de Linux.
% Y2 = log10(abs(ERB_coef_diezmados));      % Logamplitud
% fb = Y2/max(max(abs(Y2)));

% fb(fb==0) = 0.0001;              % Solución de Jesús: Sustituimos los valores de fb iguales a 0.
                                   % Con mi solución, no lo veo necesario,
                                   % puesto que no introduzco ceros.�
%--- Cálculo de los MFCCs (rastamat), de los ERBs (Auditory Toolbox + rastamat)
%    de Seneff (Auditory Toolbox + rastamat) o de dcGC (Toolbox del japonés)    ---%
ceps_ERB = fb;


%ceps_ERB = transpose(ceps_ERB);     % CPM: Jesús, he respetado aquí tu transposición porque supongo que será para algo...

%--- Cálculo de los MFCCs (Auditory Toolbox) o de los ERBs (Auditory Toolbox) ---%
%%%% ceps_ERB = log10(abs(fb));          % Logamplitud
%%%% ceps_ERB = dct(ceps_ERB); 
%%%% ceps_ERB = transpose(ceps_ERB);     % CPM: Jesús, he respetado aquí tu transposición porque supongo que será para algo...


% ceps_ERB(ceps_ERB == -Inf) = -40; % Solución de Jesús: Anulo los posibles ceps con valores -INF (claro, pues introdujo ceros).
                                    % Aunque igual funcionaría para valores de fb cercanos a 0 ¿? ¡Yo no adopto esta solución!

%--- Cálculo de los delta y delta-delta ---%
%deltaWindow    = 2;                 % Antes usaba un valor de 5 (5 tramas antes, 5 tramas después --> Erróneo).
%deltaCeps      = calculoDeltaCepstrum(ceps_ERB, deltaWindow);     % Deltas
%deltaDeltaCeps = calculoDeltaCepstrum(deltaCeps,deltaWindow);     % Doble Deltas
%c = [ceps_ERB deltaCeps deltaDeltaCeps];
c = [ceps_ERB];


end