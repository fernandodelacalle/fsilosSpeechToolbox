function y = fsilosPowSpec(x, sr, wintime, steptime, dither, SS_flag, filterbank, nbands)

%y = powspec(x, sr, wintime, steptime, sumlin, dither)
%
% compute the powerspectrum of the input signal.
% basically outputs a power spectrogram
%
% each column represents a power spectrum for a given frame
% each row represents a frequency
%
% default values:           % sr = 8000Hz
%                           wintime = 25ms (200 samps)
%                           steptime = 10ms (80 samps)
%                           which means use 256 point fft
%                           hamming window
%
% Joyner, 28/11/2011:
% 'SS_flag': 0 = off (default); 1 = on
% 'filterbank': 'mfcc', 'gammatone', 'seneff', 'dcGC'
% 'nbands': Number of frequency bands in the cochlear model.  

% Values for sr = 8000:     NFFT = 256;
%                           NOVERLAP = 120;
%                           SAMPRATE = 8000;
%                           WINDOW = hamming(200);

% Señal para verificar el flip de la matriz de coeficientes, en los 
% casos 'gammatone' o 'seneff'
% te = 0:.1:2*pi;
% x = sin(te);          % Variar la frecuencia para ver el efecto: 
                        % subida o bajada de líneas en el espectrograma

MOSTRAR = 0;            % Variable para mostrar o no las figuras.

if nargin < 2
  sr = 8000;
end
if nargin < 3
  wintime = 0.025;
end
if nargin < 4
  steptime = 0.010;
end
if nargin < 5
  dither = 1;
end
if nargin < 6
  SS_flag = 0;
end

winpts = round(wintime*sr);
steppts = round(steptime*sr);

NFFT = 2^(ceil(log(winpts)/log(2)));
WINDOW = hamming(winpts);
NOVERLAP = winpts - steppts;
SAMPRATE = sr;

if MOSTRAR
    drawspect(x,sr,1), title('Noisy Spectrogram')
end
    
% -=:=- Spectral substraction -=:=- %
% Spectral substraction seems to greatly improve the spectrogram 

if SS_flag  % SS_flag = 1 (on)      
    % Number of points in the FFT (En realidad, la mitad de los puntos
    % de la FFT: es la resoluci�n del espectrograma).
    % Valores a probar: 64, 128, 256.
    SS_resolution = 256;    
    x = fsilosSS(x, sr, SS_resolution);
    if MOSTRAR
        drawspect(x,sr,1), title('Noisy Spectrogram after Spectral Subtraction')
%%%     [Y, T2, F2] = stft(Signal_Sub,window,Overlap,fs);
%%%     figure, plotstft(Y, window, Overlap, fs), colorbar, title('Con substracci�n espectral')
    end
end

fs = SAMPRATE;                      % Sampling Frecuency
L = 512;                            % Analysis window length (equivalen a 32 mseg) 
%%% L = NFFT;          
Overlap = 352;                      % This makes the update period equal (L-Overlap) samples. (equivalen a 22 mseg)
numChannels = nbands;      
lowFreq = 100;                      % Lower frequency in the cochlear model.
%%% update_period = L - Overlap;        % (equivalen a 160 muestras, 10 mseg de desplazamiento de la ventana)

switch filterbank
    case 'mfcc'

        % Values coming out of rasta treat samples as integers, 
        % not range -1..1, hence scale up here to match (approx)
        y = abs(specgram(x*32768,NFFT,SAMPRATE,WINDOW,NOVERLAP)).^2;

    case 'gammatone'
        
        coeficientesFiltro = MakeERBFilters(fs,numChannels,lowFreq);
        %%%% ERB_coef = ERBFilterBank(x,coeficientesFiltro);
        ERB_coef = ERBFilterBank(x*32768,coeficientesFiltro);
        %%%% ERB_coef = ERBFilterBank(x*8192,coeficientesFiltro);
        
        for j = 1:size(ERB_coef,1)              % Justificaci�n de este paso en la pág. 24 de 'AuditoryToolboxTechReport.pdf'
            c = max(ERB_coef(j,:),0);           % Jesús: Coge todas las columnas de cada fila y las mete en un vector con la
                                                % condición de que los valores <0 los redondea a 0.
                                                % Joyner: Es un rectificador de media onda (valores negativos, los lleva a 0).
            %%%%    c = ERB_coef(j,:);          % Joyner: Esta línea es INCOMPATIBLE con la línea del rectificador de media onda.
            c = filter(1,[1 -0.99],c);          % Joyner: Filtro pasabajo
            ERB_coef(j,:) = c;           
        end

        %%%%% ERB_coef_diezmados = ERB_coef(:,update_period:update_period:end);
        ERB_coef_diezmados = ERBsampling(ERB_coef, WINDOW, NOVERLAP);
        ERB_coef_diezmados(ERB_coef_diezmados <= 0.001) = 0.01;         % Al principio del filtrado se ponen 0's y no desaparecen hasta llegar a
                                                                        % numChannels/2. Pongo 0.001 para evitar problemas con el logaritmo!?!?
        y = ERB_coef_diezmados;
        y = flipud(y);                                                  % Luego de hacer flipud, la frecuencia más baja queda en la primera
                                                                        % fila, y la más alta en la última fila
        if MOSTRAR
            figure, imagesc(y), title('ERBs'),colorbar                 %  Se comprueba que está bien viendo la figura después del flip
        end
        
    case 'seneff'

        x_Seneff = SeneffEar(x,fs,0);
        for j=1:40
            c=x_Seneff(j,:);
            c=filter([1],[1, -.99],c);
            x_Seneff(j,:) = c;           
        end
        Seneff_coef_diezmados = ERBsampling(x_Seneff, WINDOW, NOVERLAP);
        Seneff_coef_diezmados(Seneff_coef_diezmados <= 0.001) = 0.01;   % Al principio del filtrado se ponen 0's y no desaparecen hasta llegar a
                                                                        % numChannels/2. Pongo 0.001 para evitar problemas con el logaritmo???
        y = Seneff_coef_diezmados;                                                              
        y = flipud(y);                                                  % Luego de hacer flipud, la frecuencia más baja queda en la primera
                                                                        % fila, y la más alta en la última fila
        if MOSTRAR                                                                
            figure, imagesc(y), title('Seneffs'),colorbar              %  <--- Se comprueba que está bien viendo la figura después del flip
        end
            
    case 'dcGC'

        % Equalización del nivel de la señal
        SigSPL = 30;        % 30 dB SPL <--> RMS level of 1.0 (Meddis, 1986)
        x_eqlz =  Eqlz2MeddisHCLevel(x,SigSPL);
        %%% x_eqlz = x;     % <-- Si no se quiere equalizar la señal

        GCparam.fs     = fs;
        GCparam.NumCh  = numChannels;
        GCparam.FRange = [lowFreq, 8000];
        GCparam.OutMidCrct = 'No';      % No Outer/Middle Ear correction
        %%% GCparam.OutMidCrct = 'ELC';
        GCparam.Ctrl = 'dynamic';       % used to be 'time-varying'
        %%% GCparam.Ctrl = 'static';        % or 'fixed' 

        x_eqlz = x_eqlz';   % Es necesario un vector fila es el próx. paso
        [cGCout, pGCout, GCparam, GCresp] = GCFBv207(x_eqlz,GCparam);
        
        if MOSTRAR
            figure, imagesc(max(cGCout,0)); set(gca,'YDir','normal');
        end

        x_dcGC = ERBsampling(cGCout, WINDOW, NOVERLAP);     % cGCout:  Compressive GammaChirp Filter Output
%         x_dcGC = ERBsampling(pGCout, WINDOW, NOVERLAP);     % pGCout:  Passive GammaChirp Filter Output        
      
        x_dcGC(x_dcGC <= 0.001) = 0.01;                     % Al principio del filtrado se ponen 0's y no desaparecen hasta llegar a
                                                            % numChannels/2. Pongo 0.001 para evitar problemas con el logaritmo???
        y = x_dcGC;                                         % Con esta opción, no es necesario hacer 'flipud'.                        

        if MOSTRAR
            figure, imagesc(y), title('dcGC'),colorbar      %  Se comprueba que está bien viendo la figura después del flip        
        end
        
    otherwise
        disp(['error: Wrong filterbank value!']);
        exit;
end

% imagine we had random dither that had a variance of 1 sample 
% step and a white spectrum.  That's like (in expectation, anyway)
% adding a constant value to every bin (to avoid digital zero)
if (dither)
  y = y + winpts;
end
% ignoring the hamming window, total power would be = #pts
% I think this doesn't quite make sense, but it's what rasta/powspec.c does

% that's all she wrote
