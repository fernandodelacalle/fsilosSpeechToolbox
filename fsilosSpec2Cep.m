function [cep,dctm] = fsilosSpec2Cep(spec, ncep, type, MF_flag, se)

% [cep,dctm] = spec2cep(spec, ncep, type)
%     Calculate cepstra from spectral samples (in columns of spec)
%     Return ncep cepstral rows (defaults to 9)
%     This one does type II dct, or type I if type is specified as 1
%     dctm returns the DCT matrix that spec was multiplied by to give cep.
% 2005-04-19 dpwe@ee.columbia.edu  for mfcc_dpwe
%
% Fernando, 16/05/2013:
% 'MF_flag': 0 = off; 1 = on
% se : structural element for morfological filtering 

if nargin < 2;   ncep = 13;   end
if nargin < 3;   type = 2;   end   % type of DCT
if nargin < 4;   MF_flag = 0;   end 
if nargin < 5;   MF_flag = 0; se = 0;   end 

[nrow, ncol] = size(spec);

% Make the DCT matrix
dctm = zeros(ncep, nrow);
if type == 2 || type == 3
  % this is the orthogonal one, the one you want
  for i = 1:ncep
    dctm(i,:) = cos((i-1)*[1:2:(2*nrow-1)]/(2*nrow)*pi) * sqrt(2/nrow);
  end
  if type == 2
    % make it unitary! (but not for HTK type 3)
    dctm(1,:) = dctm(1,:)/sqrt(2);
  end
elseif type == 4 % type 1 with implicit repeating of first, last bins
  % Deep in the heart of the rasta/feacalc code, there is the logic 
  % that the first and last auditory bands extend beyond the edge of 
  % the actual spectra, and they are thus copied from their neighbors.
  % Normally, we just ignore those bands and take the 19 in the middle, 
  % but when feacalc calculates mfccs, it actually takes the cepstrum 
  % over the spectrum *including* the repeated bins at each end.
  % Here, we simulate 'repeating' the bins and an nrow+2-length 
  % spectrum by adding in extra DCT weight to the first and last
  % bins.
  for i = 1:ncep
    dctm(i,:) = cos((i-1)*[1:nrow]/(nrow+1)*pi) * 2;
    % Add in edge points at ends (includes fixup scale)
    dctm(i,1) = dctm(i,1) + 1;
    dctm(i,nrow) = dctm(i,nrow) + ((-1)^(i-1));
  end
  dctm = dctm / (2*(nrow+1));
else % dpwe type 1 - same as old spec2cep that expanded & used fft
  for i = 1:ncep
    dctm(i,:) = cos((i-1)*[0:(nrow-1)]/(nrow-1)*pi) * 2 / (2*(nrow-1));
  end
  % fixup 'non-repeated' points
  dctm(:,[1 nrow]) = dctm(:, [1 nrow])/2;
end  

%-=:=- Morphological Processing -=:=-%
if MF_flag  % MF_flag = 1 (on)      

    Y2 = log(spec);        % Logamplitud
    
    % Sin Normalización
    
    % Y2_Normalizada = Y2;
    
    % Nueva normalización (siguiendo HTK) en tres pasos):
    
    %%% Y2_Normalizada = Y2 - max(max(Y2));         % Primer paso.
    %%% Y2_Normalizada = Y2_Normalizada + 1;
    
    %%% Y2_Normalizada(Y2_Normalizada <= -4) = -4;  % Segundo paso.
    
    %%% Y2_Normalizada = 0.1*Y2_Normalizada;        % Tercer y último paso.
    
    %%% Parte de la Normalización Antigua (usada en Nolisp y en el S-Issue.
    %%%
    Y2_Normalizada = Y2/max(max(abs(Y2)));
    %%%
    
    MOSTRAR = 0;
    
    if MOSTRAR
        figure, imagesc(Y2_Normalizada), colorbar, title('Noisy Spectrogram after SS (Log-frequencies)')
        % La logfrecuencia viene dada por el uso de ERB.
    end
    
    %-=:=- Elementos Estructurantes -=:=-%
    % CPM: Elemento dise�ado para imitar el enmascaramiento del o�do humano.
    % Nhood=zeros(5,4);Nhood(:,2)=1;Nhood(2:4,3)=1;Nhood(3,4)=1;Nhood(3,1)=1;         % Structual Element (SE)
    
    %Nhood=zeros(5,7);Nhood(:,2:3)=1;Nhood(2:4,4:5)=1;Nhood(3,6:7)=1;Nhood(3,1)=1;                             % Structual Element (SE - 1er intento, utilizado para Granada y paper NoLISP)
    %se_1 = strel(Nhood);
    %%% Nhood=zeros(9,7);Nhood(5,1)=1;Nhood(:,2:3)=1;Nhood(2:8,4)=1;Nhood(3:7,5)=1;Nhood(4:6,6)=1;Nhood(5,7)=1;   % Structual Element (SE - 2do intento).
    %%% Nhood=zeros(7,7);Nhood(4,1)=1;Nhood(2:6,2)=1;Nhood(:,3:4)=1;Nhood(2:6,5)=1;Nhood(3:5,6)=1;Nhood(4,7)=1;   % Structual Element (SE - 3er intento).
    
    %Nhood=zeros(16,12);Nhood(8:9,1)=1;Nhood(7:10,2)=1;Nhood(5:12,3)=1;Nhood(3:14,4)=1;Nhood(:,5)=1;
    %Nhood(2:15,6)=1;Nhood(3:14,7)=1;Nhood(4:13,8)=1;Nhood(5:12,9)=1;Nhood(6:11,10)=1;Nhood(7:10,11)=1;Nhood(8:9,12)=1;   % Structual Element (SE - 4to intento).
    
    %%%% Nhood=zeros(7,10);Nhood(4,1)=1;Nhood(2:6,2)=1;Nhood(:,3:5)=1;Nhood(2:6,6)=1;Nhood(2:6,7)=1;Nhood(3:5,8)=1;Nhood(3:5,9)=1;Nhood(4,10)=1;   % Structual Element (SE - 5to intento).
    
    % Elementos estructurantes moviendo el centro del elemento, para que
    % coincida con el centro del enmascaramiento en tiempo y frecuencia.
    
    % Structual Element (SE - 1er intento, utilizado para Granada y paper NoLISP)
    
    % Nhood=zeros(5,9);Nhood(:,4:5)=1;Nhood(2:4,6:7)=1;Nhood(3,8:9)=1;Nhood(3,3)=1;      % Versión 1a
    
    % Nhood=zeros(5,11);Nhood(:,6:7)=1;Nhood(2:4,8:9)=1;Nhood(3,10:11)=1;Nhood(3,5)=1;   % Versión 1b
    
    % Structual Element (SE - 2do intento).
    
    %%% Nhood=zeros(9,9);Nhood(5,3)=1;Nhood(:,4:5)=1;Nhood(2:8,6)=1;Nhood(3:7,7)=1;Nhood(4:6,8)=1;Nhood(5,9)=1;   % Versión 2a
    
    %%% Nhood=zeros(9,11);Nhood(5,5)=1;Nhood(:,6:7)=1;Nhood(2:8,8)=1;Nhood(3:7,9)=1;Nhood(4:6,10)=1;Nhood(5,11)=1;  % Versión 2b
    
    % se_1 = strel(Nhood);
    
    % Elementos Estructurantes en 3D
    
    %Nhood  = zeros(5,7);Nhood(:,2:3)=1;Nhood(2:4,4:5)=1;Nhood(3,6:7)=1;Nhood(3,1)=1;        % Structual Element (SE - 1er intento, en 3D)
    %Height = zeros(5,7); Height(3,1)=0.8; Height(:,2)=0.9; Height(:,3)=1; Height(2:4,4)=1; Height(2:4,5)=0.9; Height(3,6)=0.8; Height(3,7)=0.7;
    
    %%% Nhood=zeros(9,7);Nhood(5,1)=1;Nhood(:,2:3)=1;Nhood(2:8,4)=1;Nhood(3:7,5)=1;Nhood(4:6,6)=1;Nhood(5,7)=1;   % Structual Element (SE - 2do intento).
    %%% Height = zeros(9,7); Height(5,1)=0.8; Height(:,2)=0.9; Height(:,3)=1; Height(2:8,4)=1; Height(3:7,5)=0.8; Height(4:6,6)=0.7; Height(5,7)=0.6;
    
    %%% Nhood=zeros(7,7);Nhood(4,1)=1;Nhood(2:6,2)=1;Nhood(:,3:4)=1;Nhood(2:6,5)=1;Nhood(3:5,6)=1;Nhood(4,7)=1;   % Structual Element (SE - 3er intento).
    %%% Height = zeros(7,7); Height(4,1)=0.8; Height(2:6,2)=0.9; Height(:,3)=1; Height(:,4)=1; Height(2:6,5)=0.7; Height(3:5,6)=0.6; Height(4,7)=0.5;
    
    %%% Nhood=zeros(16,12);Nhood(8:9,1)=1;Nhood(7:10,2)=1;Nhood(5:12,3)=1;Nhood(3:14,4)=1;Nhood(:,5)=1;
    %%% Nhood(2:15,6)=1;Nhood(3:14,7)=1;Nhood(4:13,8)=1;Nhood(5:12,9)=1;Nhood(6:11,10)=1;Nhood(7:10,11)=1;Nhood(8:9,12)=1;   % Structual Element (SE - 4to intento).
    %%% Height = zeros(16,12); Height(8:9,1)=0.6; Height(7:10,2)=0.7; Height(5:12,3)=0.8; Height(3:14,4)=0.9; Height(:,5)=1;
    %%% Height(2:15,6)=0.9; Height(3:14,7)=0.8; Height(4:13,8)=0.7; Height(5:12,9)=0.6; Height(6:11,10)=0.5; Height(7:10,11)=0.4; Height(8:9,12)=0.3;
    
    %%% Nhood=zeros(7,10);Nhood(4,1)=1;Nhood(2:6,2)=1;Nhood(:,3:5)=1;Nhood(2:6,6)=1;Nhood(2:6,7)=1;Nhood(3:5,8)=1;Nhood(3:5,9)=1;Nhood(4,10)=1;   % Structual Element (SE - 5to intento).
    %%% Height = zeros(7,10); Height(4,1)=0.8; Height(2:6,2)=0.9; Height(:,3:5)=1; Height(:,6)=0.9; Height(2:6,7)=0.8; Height(2:6,8)=0.7; Height(3:5,9)=0.6; Height(4,10)=0.5;
    %se_1 = strel(Nhood,Height);
    
    % La imagen utilizada para obtener la m�scara est� normalizada. Por lo
    % tanto, es "equivalente" a una imagen en escala de grises.
    
    % % Opening
    % I_erode_1 = imerode(Y2_Normalizada,se_1);
    % if MOSTRAR
    %    figure, imagesc(I_erode_1), title('Eroded 1')
    % end
    % I_dilate_1 = imdilate(I_erode_1,se_1);
    % if MOSTRAR
    %    figure, imagesc(I_dilate_1), title('Dilated 1')
    % end
    % Mask = I_dilate_1;              % M�scara logar�tmica
    
    % Closing
    I_dilate_1 = imdilate(Y2_Normalizada,se);
    if MOSTRAR
        figure, imagesc(I_dilate_1), title('Dilated 1')
    end
    I_erode_1 = imerode(I_dilate_1,se);
    if MOSTRAR
        figure, imagesc(I_erode_1), title('Eroded 1')
    end
    Mask = I_erode_1;              % M�scara logar�tmica
    
    %%% Parte de la Normalización Antigua (usada en Nolisp y en el SIssue.
    %%%
    Mask = Mask/max(max(Mask));     % Suavizado de la M�scara
    Mask(Mask <= 0.05) = 0.05;      % CPM: Limitamos en un cierto valor m�nimo
    %%%
    
    if MOSTRAR
        figure, imagesc(Mask), title('Mask')
    end
    
    %-=:=- Espectrograma Filtrado -=:=-%
    Espectrograma_Filtrado = (Y2_Normalizada + Mask)/2;  % Sumando señal ruidosa con la máscara.
    %Espectrograma_Filtrado = Mask;  % Utilizando directamente la máscara
    % para reconocimiento.
    
   % imageEspectro( Y2_Normalizada,  I_erode_1,  Espectrograma_Filtrado, SE)
    
    if MOSTRAR
        figure, imagesc(Espectrograma_Filtrado), title('Noisy Spectrogram after SS and MF'),colorbar
    end
    
    cep = dctm*Espectrograma_Filtrado;

else % MF_flag = 0 (off)
    cep = dctm*log(spec);   

end