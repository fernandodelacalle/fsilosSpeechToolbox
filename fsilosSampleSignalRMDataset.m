function [s , fs] = fsilosSampleSignalRMDataset()
    audio = 'SampleFiles/sampleFromRMDataset.wav';
    [ s , fs ] = fsilosReadRM( audio );
end