function [SE] = se3dHyperboloid(t_fin, t_ini, time_step, N)
% function [SE,M_temp,M_frec] = se3dHyperboloid(t_fin, t_ini, time_step, frequencySize)
% 
%Returns

    error(nargchk(4, 4, nargin));
    
    %% Create the two zones
    DeltaTBack = t_ini:time_step:0;
    DeltaTFor = 0:time_step:t_fin;
    
    centralBand = 0;
    DeltaF = centralBand:N;
    DeltaF = [centralBand-N:centralBand, DeltaF(2:end)];
    %[~, DeltaF] = frecuencyMaskingNew(frequencySize, 0);
    %[M_temp , DeltaT] = temporalMaskingNew(t_fin,t_ini, time_step);
    
    %Define the hyperboloid: Lower sheet is selected by the minus sign
    hypp = @(x,y,a,b,c)  - (c .* sqrt(1 + x.^2./a^2 + y.^2./b^2));
    
    %% Define four regions
    %% Forward amount of masking
    [x,y] = meshgrid(DeltaTFor,DeltaF);
    THQFor = zeros(length(DeltaF),length(DeltaTFor));%Nice way to do it, mainly
    %readjust the forward paraboloid
    c = 5;
    SEFor = max(THQFor,4*c + hypp(x,y,50,3,c));
    
    %% Backward amount of masking
    [x,y] = meshgrid(DeltaTBack,DeltaF);
    THQBack = zeros(length(DeltaF),length(DeltaTBack));%Nice way to do it, mainly
    SEBack = max(THQBack,4*c + hypp(x,y,50,3,c));
    
    %PAdding 
    addZeros = ((t_fin + t_ini) / time_step) - 1;
    SE = [zeros(N*2 + 1,addZeros),  SEBack, SEFor];
    %Normalize
    SE = (SE - min(min(SE))) ./ (max(max(SE)) - min(min(SE)));
end
