% [ad_y] = fsilosRecoverFromPNCC(ad_x, fs, SS_flag)
% Programmed by Fernando de la Calle based in the Chanwoo Kim
% implementation of the PNCC
% Important: The input should be mono and be sampled at "16 kHz".
function [ad_y] = fsilosRecoverFromPNCC(ad_x, fs, SS_flag)
    if nargin < 2;   error('s and fs requiered'); return;  end      
    if nargin < 3;   SS_flag = 0;   end    
    % Spectral substraction
    % Spectral substraction seems to greatly improve the spectrogram 
    if SS_flag  % SS_flag = 1 (on)      
        % Number of points in the FFT (En realidad, la mitad de los puntos
        % de la FFT: es la resolucion del espectrograma).
        % Valores a probar: 64, 128, 256.
        SS_resolution = 256;    
        ad_x = fsilosSS(ad_x, fs, SS_resolution);
    end   
    % Flags
	bPreem         = 1;
    % Constants
	dSampRate   = fs;
	dLowFreq      = 200;
	dHighFreq     = dSampRate / 2;
	dFactor = 2.0;
	dGammaThreshold = 0.005;
	iM = 2;
	iN = 4;
	iSMType = 0;    
	dLamda  = 0.999;
	dLamda2 = 0.5;
	dDelta1 = 1;
	dLamda3 = 0.85;
	dDelta2 = 0.2;
	dFrameLen     = 0.0256;  % 25.6 ms window length, which is the default setting in CMU Sphinx
	dFramePeriod = 0.010;   % 10 ms frame period
	iFFTSize  = 1024;
	iNumFilts = 40;
	% Array Queue Ring-buffer
	global Queue_aad_P;
	global Queue_iHead;
	global Queue_iTail;
	global Queue_iWindow;
	global Queue_iNumElem;
	Queue_iWindow  = 2 * iM + 1;
	Queue_aad_P    = zeros(Queue_iWindow, iNumFilts);
	Queue_iHead    = 0;
	Queue_iTail    = 0;
	Queue_iNumElem = 0;
    %
	iFL        = floor(dFrameLen    * dSampRate);
	iFP        = floor(dFramePeriod * dSampRate);
	iNumFrames = floor((length(ad_x) - iFL) / iFP) + 1;
	iSpeechLen = length(ad_x);
	% Pre-emphasis using H(z) = 1 - 0.97 z ^ -1
    if (bPreem == 1)
	    ad_x = filter([1 -0.97], 1, ad_x);
    end
	% Obtaning the gammatone coefficient. 
	% Based on M. Snelly's auditory toolbox. 
	% In actual C-implementation, we just use a table
   	aad_H = fsilosComputeFilterResponse(iNumFilts, iFFTSize, dLowFreq, dHighFreq, dSampRate);
	aad_H = abs(fsilosNormalizeFilterGain(aad_H, dSampRate));
	for i = 1 : iNumFilts,
	    aiIndex = find(aad_H(:, i) < dGammaThreshold * abs(max(aad_H(:, i))));
	    aad_H(aiIndex, i) = 0;
    end 
    % Preallocation variables
	aad_P           = zeros(iNumFrames,      iNumFilts);
	ad_Q            = zeros(1,               iNumFilts);
	ad_Q_Out        = zeros(1,               iNumFilts);
	ad_QMVAvg       = zeros(1,               iNumFilts);
	ad_w            = zeros(1,               iNumFilts);
	ad_w_sm         = zeros(1,               iNumFilts);
	ad_QMVAvg_LA    = zeros(1,               iNumFilts);     
    ad_y            = zeros(size(ad_x));
	% Obtaining the short-time Power P(i, j)
    i_FI     = 0;
    for m = 0 : iFP : iSpeechLen  - iFL 
        ad_x_st                = ad_x(m + 1 : m + iFL) .* hamming(iFL);
        adSpec                 = fft(ad_x_st, iFFTSize);
        ad_X                   = abs(adSpec(1: iFFTSize / 2));
        aadX(:, i_FI + 1)      = ad_X;        
        ad_mu                  = zeros(1, iFFTSize / 2);
        % Calculating the Power P(i, j)
        for j = 1 : iNumFilts
            % Squared integration
            aad_P(i_FI + 1, j)  = sum((ad_X .* aad_H(:, j)) .^ 2);
            	        
        end
        % Ring buffer (using a Queue)
        if (i_FI >= 2 * iM + 1)
            Queue_poll();
        end
        Queue_offer(aad_P(i_FI + 1, :));
        ad_Q = Queue_avg();
        if (i_FI == 2 * iM)
            ad_QMVAvg     = ad_Q.^ (1 / 15);
            ad_PBias  =  (ad_Q) * 0.9;
        end         
        if (i_FI >= 2 * iM)  
            % Bias Update
            for i = 1 : iNumFilts,
                if (ad_Q(i) > ad_PBias(i))
                   ad_PBias(i) = dLamda * ad_PBias(i)  + (1 - dLamda) * ad_Q(i);
                else
                   ad_PBias(i) = dLamda2 * ad_PBias(i) + (1 - dLamda2) * ad_Q(i);
                end
            end		            
            for i = 1 : iNumFilts,
                ad_Q_Out(i) =   max(ad_Q(i) - ad_PBias(i), 0) ;
                if (i_FI == 2 * iM)
                    ad_QMVAvg2(i)  =  0.9 * ad_Q_Out(i);
                    ad_QMVAvg3(i)  =  ad_Q_Out(i);
                    ad_QMVPeak(i)  =  ad_Q_Out(i);
                end
                if (ad_Q_Out(i) > ad_QMVAvg2(i))
                     ad_QMVAvg2(i) = dLamda * ad_QMVAvg2(i)  + (1 -  dLamda)  *  ad_Q_Out(i);
                else
                     ad_QMVAvg2(i) = dLamda2 * ad_QMVAvg2(i) + (1 -  dLamda2) *  ad_Q_Out(i);
                end
                dOrg =  ad_Q_Out(i);
                ad_QMVAvg3(i) = dLamda3 * ad_QMVAvg3(i);
                if (ad_Q(i) <  dFactor * ad_PBias(i))
                    ad_Q_Out(i) = ad_QMVAvg2(i);
                else
                     if (ad_Q_Out(i) <= dDelta1 *  ad_QMVAvg3(i))
                        ad_Q_Out(i) = dDelta2 * ad_QMVAvg3(i);
                     end
                end
                ad_QMVAvg3(i) = max(ad_QMVAvg3(i),   dOrg);
                ad_Q_Out(i) =  max(ad_Q_Out(i), ad_QMVAvg2(i));
            end
            ad_w      =   ad_Q_Out ./ max(ad_Q, eps);
            for i = 1 : iNumFilts,
                if iSMType == 0
                        ad_w_sm(i) = mean(ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts)));
                elseif iSMType == 1
                        ad_w_sm(i) = exp(mean(log(ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts)))));
                elseif iSMType == 2
                        ad_w_sm(i) = mean((ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts))).^(1/15))^15;
                elseif iSMType == 3
                        ad_w_sm(i) = (mean(  (ad_w(max(i - iN, 1) : min(i + iN ,iNumFilts))).^15 )) ^ (1 / 15); 
                end
            end
            for channelsaux = 1 : 40
                ad_mu = ad_mu + ad_w_sm(channelsaux) * abs(aad_H(:, channelsaux))';        
            end
            % Denominator part of Equation (7) in the IS2010 paper
            ad_mu      = ad_mu ./ (sum(abs(aad_H')  + eps)); 
            % Equation (8) in the IS2010 paper
            ad_mu_sym  = [ad_mu, fliplr(ad_mu)];  
            % Equation (9) in the IS2010 paper
            ad_y_st               = ifft(adSpec .* ad_mu_sym');
            % Overlap addition to resynthesize speech
            ad_y(m + 1 : m + iFL) = ad_y(m + 1 : m + iFL) + ad_y_st(1 : iFL);        
        end
        i_FI = i_FI + 1;
    end
    % Prepare the final signal
    ad_y = real(ad_y);
    ad_y = filter(1, [1 -0.97], ad_y);
    % Scaling and writing the file
    dStdOrg         = std(ad_x);
    ad_y = (ad_y - mean(ad_y)) * dStdOrg / std(ad_y);
    dMaxAmp = max(abs(ad_y));
    % To prevent clipping
    if (dMaxAmp >= 1)
        ad_y = ad_y /dMaxAmp / 1.01; 
    end
end

function [] = Queue_offer(ad_x)
    global Queue_aad_P;
    global Queue_iHead;
    global Queue_iTail;
    global Queue_iWindow;
    global Queue_iNumElem;
    Queue_aad_P(Queue_iTail + 1, :) = ad_x;    
    Queue_iTail    = mod(Queue_iTail + 1, Queue_iWindow);
    Queue_iNumElem = Queue_iNumElem + 1;
    if Queue_iNumElem > Queue_iWindow
       error ('Queue overflow'); 
    end  
end

function [ad_x] = Queue_poll()
    global Queue_aad_P;
    global Queue_iHead;
    global Queue_iTail;
    global Queue_iWindow;
    global Queue_iNumElem;
    if Queue_iNumElem <= 0
       error ('No elements'); 
    end
    ad_x =  Queue_aad_P(Queue_iHead + 1, :);   
    Queue_iHead    = mod(Queue_iHead + 1, Queue_iWindow);
    Queue_iNumElem = Queue_iNumElem - 1;
end

function[adMean] = Queue_avg()
    global Queue_aad_P;
    global Queue_iHead;
    global Queue_iTail;
    global Queue_iWindow;
    global Queue_iNumElem;  
    adMean = zeros(1, 40);
    iPos = Queue_iHead;  
    for i = 1 : Queue_iNumElem
        adMean = adMean + Queue_aad_P(iPos + 1 ,: );
        iPos   = mod(iPos + 1, Queue_iWindow);
    end   
    adMean = adMean / Queue_iNumElem;
end
