function [ss_noisy] =  fsilosAddNoiseFANT(ss,tipoRuido, snr)        
    auxVariables =fsilosVariables();
    dirfsilosToolbox = auxVariables{1};
    ruidos = {'OriginalNoisesTransPaperPNCC/whitenoise.raw','OriginalNoisesTransPaperPNCC/streetNoise.raw','OriginalNoisesTransPaperPNCC/music.raw','OriginalNoisesTransPaperPNCC/interfering.raw'};
    rawFile = 'fsilosTemp.wav'; 
    fsilosWriteRawFile(ss, rawFile, 'int16');
    fich_in = fopen('fsilosFATNin.list','wt');
    fprintf(fich_in,rawFile);
    fclose(fich_in);    
    ruidoDir = [dirfsilosToolbox  ruidos{tipoRuido}];
    filterAddNoiseCommand = [dirfsilosToolbox 'filter_add_noise'];    
    system([filterAddNoiseCommand ' -u -i fsilosFATNin.list -o fsilosFATNin.list -n ' ruidoDir ' -s ' num2str(snr)]);
    ss_noisy = fsilosReadRawfile(rawFile,'int16');
    delete('fsilosFATNin.list');
    delete(rawFile);
    delete('filter_add_noise.log')  
end