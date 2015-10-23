%
% Programmed by  Chanwoo Kim
%
% Read feature in sphinx format.
%
% The first 4-byte is in 
%
% The default value of iDim is 13...
function [aadFeat] = fsilosReadSphinxFeat(szFileName, iDim)
% Sample feat file... 
%
%E:\RM_WSJ_NO_PROCESSING_BASE_0_4_ICASSP\RM_test1600_White_10dB\aem0_2\sb02
%.mfc
%

if nargin < 2
    iDim = 13;
end

fid = fopen(szFileName, 'rb');

if (fid <= 0)
    error('file cannot be opened.');
end

bHeader = fread(fid, 1, 'int32');
aadFeat = fread(fid, [iDim, inf], 'float32');
[M, N] = size(aadFeat);
bHeader
        M 
        N
if (bHeader ~= M * N)
    fseek(fid, 0, 'bof');
    bHeader = fread(fid, 1, 'int32', 'ieee-be');
    aadFeat = fread(fid, [iDim, inf], 'float32', 'ieee-be');
    [M, N] = size(aadFeat);    
    if (bHeader ~= M * N)
        bHeader
        M 
        N
        error ('size does not match') 
    end 
end

end
    


