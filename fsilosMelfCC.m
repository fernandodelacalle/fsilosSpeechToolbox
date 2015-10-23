function [cepstra,aspectrum,pspectrum] = fsilosMelfCC(samples, sr, filterbank, MF_flag, se,  SS_flag, varargin) 

%[cepstra,aspectrum,pspectrum] = melfcc(samples, sr[, opts ...])
%  Calculate Mel-frequency cepstral coefficients by:
%   - take the absolute value of the STFT
%   - warp to a Mel frequency scale
%   - take the DCT of the log-Mel-spectrum
%   - return the first <ncep> components
%  This version allows a lot of options to be controlled, as optional 
%  'name', value pairs from the 3rd argument on: (defaults in parens)
%    'wintime' (0.025): window length in sec
%    'hoptime' (0.010): step between successive windows in sec
%    'numcep'     (13): number of cepstra to return
%    'lifterexp' (0.6): exponent for liftering; 0 = none; < 0 = HTK sin lifter
%    'sumpower'    (1): 1 = sum abs(fft)^2; 0 = sum abs(fft)
%    'preemph'  (0.97): apply pre-emphasis filter [1 -preemph] (0 = none)
%    'dither'      (0): 1 = add offset to spectrum as if dither noise
%    'minfreq'     (0): lowest band edge of mel filters (Hz)
%    'maxfreq'  (4000): highest band edge of mel filters (Hz)
%    'nbands'     (40): number of warped spectral bands to use
%    'bwidth'    (1.0): width of aud spec filters relative to default
%    'dcttype'     (2): type of DCT used - 1 or 2 (or 3 for HTK or 4 for feac)
%    'fbtype'  ('mel'): frequency warp: 'mel','bark','htkmel','fcmel'
%    'usecmp'      (0): apply equal-loudness weighting and cube-root compr.
%    'modelorder'  (0): if > 0, fit a PLP model of this order
% The following non-default values nearly duplicate Malcolm Slaney's mfcc
% (i.e. melfcc(d,16000,opts...) =~= log(10)*2*mfcc(d*(2^17),16000) )
%       'wintime': 0.016
%     'lifterexp': 0
%       'minfreq': 133.33
%       'maxfreq': 6855.6
%      'sumpower': 0
% The following non-default values nearly duplicate HTK's MFCC
% (i.e. melfcc(d,16000,opts...) =~= 2*htkmelfcc(:,[13,[1:12]])'
%  where HTK config has PREEMCOEF = 0.97, NUMCHANS = 20, CEPLIFTER = 22, 
%  NUMCEPS = 12, WINDOWSIZE = 250000.0, USEHAMMING = T, TARGETKIND = MFCC_0)
%     'lifterexp': -22
%        'nbands': 20
%       'maxfreq': 8000
%      'sumpower': 0
%        'fbtype': 'htkmel'
%       'dcttype': 3
% For more detail on reproducing other programs' outputs, see
% http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/mfccs.html
%
% 2005-04-19 dpwe@ee.columbia.edu after rastaplp.m.  
% Uses Mark Paskin's process_options.m from KPMtools
%
% Joyner, 28/11/2011:
% 'filterbank': 'mfcc', 'gammatone', 'seneff', 'dcGC'
% 'MF_flag': 0 = off, 1 = on, Morphological Filtering.

if nargin < 2;   sr = 16000;    end

% Parse out the optional arguments
[wintime, hoptime, numcep, lifterexp, sumpower, preemph, dither, ...
 minfreq, maxfreq, nbands, bwidth, dcttype, fbtype, usecmp, modelorder] = ...
    process_options(varargin, 'wintime', 0.025, 'hoptime', 0.010, ...
          'numcep', 13, 'lifterexp', 0.6, 'sumpower', 1, 'preemph', 0.97, ...
	  'dither', 0, 'minfreq', 0, 'maxfreq', 4000, ...
	  'nbands', 40, 'bwidth', 1.0, 'dcttype', 2, ...
	  'fbtype', 'mel', 'usecmp', 0, 'modelorder', 0);

if preemph ~= 0
  samples = filter([1 -preemph], 1, samples);
end

% Compute FFT power spectrum
pspectrum = fsilosPowSpec(samples, sr, wintime, hoptime, dither, SS_flag, filterbank, nbands);
           
switch filterbank
    case 'mfcc'
        aspectrum = audspec(pspectrum, sr, nbands, fbtype, minfreq, maxfreq, sumpower, bwidth);        

        % Problema detectado: Con 128 bandas, la primera fila del aspectrum
        % se convierte en una fila de ceros (corresponde a la frecuencia
        % más baja). Luego, al calcular el log en el script "scep2cep" 
        % aparecen "-Inf". A continuación sustituimos esa fila por un valor 
        % muy pequeño cercano a 0 (así el log será un número negativo 
        % tendiendo a -infinito. Otra alternativa es sustituir el resultado
        % del log, -Inf, por un número negativo grande (cercano a -Inf).
        
        % Es fácil darse cuenta que es la frecuencia más baja mirando en el
        % script powspec.m el switch de filterbank, y los flipud de otros
        % casos.
        
        if nbands == 128
            aspectrum(1,:) = 0.000000000000000001;
        end
        
    case {'gammatone', 'seneff', 'dcGC'}
        aspectrum = pspectrum;      % En estos casos, no hace falta hacer aparte el análisis de bandas        
    
    otherwise
        disp(['error: Wrong filterbank value!']);
        exit;
end
   
if (usecmp)
  % PLP-like weighting/compression
  aspectrum = postaud(aspectrum, maxfreq, fbtype);
end

if modelorder > 0

  if (dcttype ~= 1) 
    disp(['warning: plp cepstra are implicitly dcttype 1 (not ', num2str(dcttype), ')']);
  end
  
  % LPC analysis 
  lpcas = dolpc(aspectrum, modelorder);

  % convert lpc to cepstra
  cepstra = lpc2cep(lpcas, numcep);

  % Return the auditory spectrum corresponding to the cepstra?
%  aspectrum = lpc2spec(lpcas, nbands);
  % else return the aspectrum that the cepstra are based on, prior to PLP

else
  
  % Convert to cepstra via DCT
  cepstra = fsilosSpec2Cep(aspectrum, numcep, dcttype, MF_flag, se);

end

cepstra = lifter(cepstra, lifterexp);

