function fsilosCompareTwoSignals(s1, s2, fs)
    %ss = ss ./ max(ss);
    soundsc(s1,fs);
    %ss_noisy = ss_noisy ./ max(ss_noisy);
    soundsc(s2,fs)   
    
    figure(1);
    subplot(211);
    plot(s1);
    subplot(212);
    plot(s2);     
end
