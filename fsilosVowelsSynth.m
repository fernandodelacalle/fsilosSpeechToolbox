% [xout , fs] = fsilosVowelsSynth(s, fs, f0, dur)
% No parameters signal will be recorded from mic
% f0 in Hz, excitation frequency 
% dur = 3 duration in sec
% Simple script to record vowels and synthesize in monotone form using LPC
% written 10/13/15 - rms
function [xout , fs] = fsilosVowelsSynth(s, fs, f0, dur, playsound)
    record = 0;
    if nargin < 2;   
        record = 1;  
    end   
    if nargin < 3;       
        f0 = 100;  % in Hz, excitation frequency 
    end
    if nargin < 4;   
        dur = 1;  
    end 
    if nargin < 5;
        playsound = 0;
    end
    
 
    if record
        fs = 16000;   
        % generate excitation
        excite = zeros(fs*dur,1);
        for n = 1:length(excite)
            if mod(n,round(fs/f0)) == 0
                excite(n) = 1;
            end
        end
        % Record vowel
        n = 0;  
        while n == 0
            recobj = audiorecorder(fs,16,1);
            disp('Hit return to begin')
            pause
            recordblocking(recobj,dur);
            disp('Recording over')
            play(recobj)
            reply = input('Recording OK [y/n]?','s')
            if reply == 'y'
                n = 1;
            end
        end
        xraw = getaudiodata(recobj);
    else  
        % generate excitation
        excite = zeros(fs*dur,1);
        for n = 1:length(excite)
            if mod(n,round(fs/f0)) == 0
                excite(n) = 1;
            end
        end
        xraw = s;
    end
    
    a = lpc(xraw,14);
    xout = filter(1,a,excite);
    if playsound
        soundsc(xout,fs);
    end
    freqz(1,a);
end