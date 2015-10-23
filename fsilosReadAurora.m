function [s, fs] = fsilosReadAurora(name)
    fs = 8000;
    fid = fopen(name, 'r','ieee-be');
    s = fread(fid, 'int16');
    % Normalizacion agregada el 5/11/2013.
    s = s/max(abs(s)); 
    fclose(fid);
end