%
% Programmed by Chanwoo Kim for the IEEETran Speech, Audio, and Langauge Processing and ICASSP 2012
%
% (chanwook@cs.cmu.edu)
%
% Important: The input should be mono and be sampled at "16 kHz".
%
%
% * If you want to use logarithmic nonlinearity instead of the power
% nonlinearity, change bPowerLaw to 0 (lilne 28)
%
% PNCC_IEEETran(OutFile, InFile)
%
% iFiltTyep 2 Gamma
% iFiltType 1 HTK MEL
% iFiltType 0 Snaley MEL
% Default
% 0.5, 0.01, 2, 4
%
function [aadDCT] = PNCC__sincFer_v_0_1(s, fs, MF_flag, se)

% 	fid = fopen(szInFileName, 'rb');
% 	fseek(fid, 1024, 'bof');
% 	ad_x  = fread(fid, 'int16');
% 	fclose(fid);

    ad_x = s;
    
	dLamda_L = 0.999;
	dLamda_S = 0.999;

	dSampRate   = fs;
	dLowFreq      = 200;
	dHighFreq     = dSampRate / 2;
	dPowerCoeff = 1 / 15;

	iFiltType = 2;
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
  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Flags
    %
	bPreem         = 0;
	bSSF             = 1;
	bPowerLaw    = 0;
	bDisplay        = 0;
     

	dFrameLen     = 0.0256;  % 25.6 ms window length, which is the default setting in CMU Sphinx
	dFramePeriod = 0.010;   % 10 ms frame period
	iPowerFactor  = 1;

	iFFTSize  = 1024;
	iNumFilts = 40;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Flags
	%
	%
	% Array Queue Ring-buffer
	%
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
   
	iFL        = floor(dFrameLen    * dSampRate);
	iFP        = floor(dFramePeriod * dSampRate);
	iNumFrames = floor((length(ad_x) - iFL) / iFP) + 1;
	iSpeechLen = length(ad_x);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Pre-emphasis using H(z) = 1 - 0.97 z ^ -1
	%
	if (bPreem == 1)
	    ad_x = filter([1 -0.97], 1, ad_x);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Obtaning the gammatone coefficient. 
	%
	% Based on M. Snelly's auditory toolbox. 
	% In actual C-implementation, we just use a table
	%
	
    
	i_FI     = 0;
	i_FI_Out = 0;

	if bSSF == 1
	    adSumPower = zeros(1, iNumFrames - 2 * iM);
	else
	    adSumPower = zeros(1, iNumFrames);
	end
	 
	%dLamda_L   = 0.998;
	aad_P      = zeros(iNumFrames,      iNumFilts);
	aad_P_Out  = zeros(iNumFrames - 2 * iM,      iNumFilts);
	ad_Q       = zeros(1,               iNumFilts);
	ad_Q_Out   = zeros(1,               iNumFilts);
	ad_QMVAvg  = zeros(1,               iNumFilts);
	ad_w       = zeros(1,               iNumFilts);
	ad_w_sm    = zeros(1,               iNumFilts);
	ad_QMVAvg_LA = zeros(1,               iNumFilts);

	MEAN_POWER = 1e10;

	dMean  = 5.8471e+08;
	dPeak = 2.7873e+09 / 15.6250;
	% (1.7839e8, 2.0517e8, 2.4120e8, 2.9715e8, 3.9795e8) 95, 96, 97, 98, 99
	% percentile from WSJ-si84
	            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dPeakVal = 4e+07;% % 4.0638e+07  --> Mean from WSJ0-si84  (Important!!!)
	                %%%%%%%%%%%
	dMean = dPeakVal;
    
	
    
    % Sennef Model
    bm = SeneffEar(ad_x, fs, 0);
    % Sincrony detector
    [bm_ret] = (GSD(bm, fs));  
    bm_ret = bm_ret +1e-2;           
    % LPF
    i_FI = 0;
    for m = 0 : iFP : iSpeechLen  - iFL 
        aad_P(i_FI + 1,:) = sum(bm_ret(:,(m + 1 : m + iFL),:),2);
        adSumPower(i_FI + 1) = sum(aad_P(  i_FI + 1, :));
        i_FI = i_FI + 1;
    end
    
    % Peak Power Normalization Using 95 % percentile
% 	adSorted  = sort(adSumPower);
% 	dMaxPower = adSorted(round(0.95 * length(adSumPower)));
%	aad_P     = aad_P / dMaxPower * 1e15;
     
   ch = 8;
   figure();
   subplot(411);
   plot(bm(ch,:));
   subplot(412);
   plot(bm_ret(ch,:));
   subplot(413);   
   plot(aad_P(:,ch));
   subplot(414);   
   plot(adSumPower);


    
    
    i_FI = 0;
    for m = 0 : iFP : iSpeechLen  - iFL 
        if bSSF == 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Ring buffer (using a Queue)
                %
                if (i_FI >= 2 * iM + 1)
                    Queue_poll();
                end
                Queue_offer(aad_P(i_FI + 1, :));
                ad_Q = Queue_avg();

                %FER
                Q_avg(i_FI + 1 ,:) = ad_Q;
                
                if (i_FI == 2 * iM)
                    ad_QMVAvg     = ad_Q.^ (1 / 15);
                    ad_PBias  =  (ad_Q) * 0.99999;
                end

                if (i_FI >= 2 * iM)  
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % Bias Update
                    %
                    for i = 1 : iNumFilts,
                        if (ad_Q(i) > ad_PBias(i))
                           ad_PBias(i) = dLamda * ad_PBias(i)  + (1 - dLamda) * ad_Q(i);
                        else
                           ad_PBias(i) = dLamda2 * ad_PBias(i) + (1 - dLamda2) * ad_Q(i);
                        end      
                    end
                    
                    % FER
                    Q_bias(i_FI - iM + 1, :) = ad_PBias(:);
                    
                    
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
                        
                        %FER
                        Q_out(i_FI - iM + 1, i) = ad_Q_Out(i);
                        Q_in(i_FI - iM + 1 ,i) = ad_QMVAvg2(i); 
                        
                        
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

                    aad_P_Out(i_FI_Out + 1, :) = ad_w_sm .* aad_P(i_FI - iM + 1, :);
    % 		        adSumPower(i_FI_Out + 1)   = sum(aad_P_Out(i_FI_Out + 1, :));
    % 
    % 		        if  adSumPower(i_FI_Out + 1) > dMean
    % 		             dMean = dLamda_S * dMean + (1 - dLamda_S) * adSumPower(i_FI_Out + 1);
    % 		        else
    % 		             dMean = dLamda_L * dMean + (1 - dLamda_L) * adSumPower(i_FI_Out + 1);
    % 		        end
    % 		        
    % 		        aad_P_Out(i_FI_Out + 1, :) = aad_P_Out(i_FI_Out + 1, :) / (dMean)  * MEAN_POWER;
                    i_FI_Out = i_FI_Out + 1;		        
                end
           
        else % if not SSF		
%         % Mean Power normalization
%         adSumPower(i_FI + 1)   = sum(aad_P(i_FI + 1, :));     
%         % dLamda_S = dLamda_L = 0.999
%         if  adSumPower(i_FI_Out + 1) > dMean
% 		     dMean = dLamda_S * dMean + (1 - dLamda_S) * adSumPower(i_FI_Out + 1);
% 		  else
% 		     dMean = dLamda_L * dMean + (1 - dLamda_L) * adSumPower(i_FI_Out + 1);
% 		  end
            aad_P_Out(i_FI + 1, :) = aad_P(i_FI + 1, :);
        end
        i_FI = i_FI + 1;
    end

    
    
    ch = 8;
    figure();
    subplot(411);
    plot(aad_P(:, ch));
    subplot(412);
    plot(Q_avg(:, ch));   
    subplot(413);
    plot(Q_bias(:,ch));
    subplot(414);   
    plot(Q_out(:,ch));

    
    figure();
    hold all
    plot(Q_avg(5:end, ch));   
    plot(Q_bias(5:end,ch));
    
    
    % Peak power normalization
    %adSorted  = sort(adSumPower);
    %dMaxPower = adSorted(round(0.95 * length(adSumPower)));
    %aad_P_Out     = aad_P / dMaxPower * 1e15;
	% Apply the nonlinearity
    if bPowerLaw == 1
	    aadSpec = aad_P_Out .^ dPowerCoeff;
	else
	    %aadSpec = log(aad_P_Out + eps);
        aadSpec = (aad_P);
    end
	% MORPHOLOGICAL FILTER
    if MF_flag
        Y2_Normalizada = aadSpec;
        I_dilate_1 = imdilate(Y2_Normalizada,se);
        I_erode_1 = imerode(I_dilate_1,se);
        Mask = I_erode_1;              
        %%% Parte de la Normalización Antigua (usada en Nolisp y en el SIssue.
        Mask = Mask/max(max(Mask));     % Suavizado de la M�scara
        Mask(Mask <= 0.05) = 0.05;      % CPM: Limitamos en un cierto valor m�nimo
        Espectrograma_Filtrado = (Y2_Normalizada + Mask)/2;  % Sumando señal ruidosa con la máscara.   
    else
        Espectrograma_Filtrado = aadSpec;
    end
	% DCT
	aadDCT                  = dct(Espectrograma_Filtrado')';
	aadDCT(:, 14:iNumFilts) = [];
	% CMN
	for i = 1 : 13
	       aadDCT( :, i ) = aadDCT( : , i ) - mean(aadDCT( : , i));
    end
	aadDCT = aadDCT';
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

function [sincrony] = GSD(synapse, fs)  
    


FilterBankRTheta = [
     0         3.14159       0.740055     2.633909    0.8
     0.86      2.997077      0.753637     2.178169    0.8
     0.86      2.879267      0.775569     1.856744    0.8
     0.86      2.761458      0.798336     1.617919    0.8
     0.86      2.643648      0.819169     1.433496    0.8
     0.86      2.525839      0.837158     1.286795    0.8
     0.8       2.964876      0.852598     1.167321    0.8
     0.86      2.408029      0.865429     1.068141    0.8
     0.86      2.29022       0.876208     0.984489    0.8
     0.86      2.17241       0.885329     0.912985    0.8
     0.86      2.054601      0.893116     0.851162    0.8
     0.86      1.936791      0.899823     0.797179    0.8
     0.8       2.788161      0.906118     0.749633    0.8
     0.86      1.818981      0.911236     0.70744     0.8
     0.86      1.701172      0.915747     0.669742    0.8
     0.86      1.583362      0.919753     0.635858    0.8
     0.86      1.465552      0.923335     0.605237    0.8
     0.86      1.347743      0.926565     0.57743     0.8
     0.8       2.611447      0.929914     0.552065    0.8
     0.86      1.229933      0.932576     0.528834    0.8
     0.86      1.112123      0.944589     0.487783    0.75
     0.86      0.994314      0.957206     0.452645    0.660714
     0.86      0.876504      0.956548     0.42223     0.672143
     0.86      0.758694      0.956653     0.395644    0.682143
     0.8       2.434732      0.956518     0.372208    0.690966
     0.86      0.640885      0.956676     0.351393    0.69881
     0.86      0.523075      0.956741     0.316044    0.712143
     0.8       2.258018      0.956481     0.287157    0.723052
     0.8       2.081304      0.956445     0.263108    0.732143
     0.8       1.904589      0.956481     0.242776    0.739835
     0.86      0.405265      0.958259     0.217558    0.749384
     0.8       1.727875      0.963083     0.197086    0.757143
     0.8       1.55116       0.969757     0.175115    0.769048
     0.8       1.374446      0.97003      0.153697    0.780662
     0.8       1.197732      0.970382     0.134026    0.791337
     0.8       1.021017      0.970721     0.118819    0.799596
     0.8       1.5           0.970985     0.106711    0.8
     0.8       1.2           0.971222     0.096843    0.8
     0.8       1             0.97144      0.088645    0.8
     0.8       0.9           0.971645     0.081727    0.8];

 
    cfs = FilterBankRTheta(:,4)/pi*fs/2;
    n_i = round(fs ./cfs);
    beta = 0.9999;
    delta = 0.01;
    As = 3;
   
    a = zeros(size(synapse));
    b = zeros(size(synapse));
    for i = 1:40
       delayFilter = dfilt.delay(n_i(i)); 

       celloutput = synapse(i, :);
       delaySignal = filter(delayFilter, celloutput);

       a(i,:)  = abs(celloutput + delaySignal);
       b(i,:)  = abs(celloutput - beta*delaySignal);
    end
    
    
    aa = 50/fs;
    a_mr = abs(filter(aa, [1 aa-1], a, [], 2));
    b_mr = abs(filter(aa, [1 aa-1], b, [], 2));

    a_mr = ( a_mr - delta);
    
    sincrony = As .* atan( (1/As)*(a_mr ./b_mr)  );
    
    a_mr = a_mr - delta;
    
    sincrony = As .* atan( (1/As)*(a_mr ./ b_mr)  );
end


