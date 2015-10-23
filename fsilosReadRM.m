function [ s , fs ] = fsilosReadRM( audioFile )
    fs = 16000;   
    fid = fopen(audioFile, 'rb');
    fseek(fid, 1024, 'bof');
    s  = fread(fid, 'int16');
    fclose(fid);    
end

