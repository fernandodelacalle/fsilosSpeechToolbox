% function fsilosWriteSphinxFeature( c , OutFeatFileName)
function fsilosWriteSphinxFeature( c , OutFeatFileName)
    [iM, iN] = size(c);
    iNumData = iM * iN;
    fid = fopen(OutFeatFileName, 'wb');
    fwrite(fid, iNumData, 'int32');
    iCount = fwrite(fid, c(:), 'float32');
    fclose(fid);
end