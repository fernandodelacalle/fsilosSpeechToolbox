function [SE] = se3dAsymHyperboloid(t_fin, t_ini, time_step, N,T1,T2,F1,F2)
% function [SE] = se3dAsymHyperboloid(t_fin, t_ini, time_step, N,T1,T2,F2,F2)
%
% Returns a structuring element with approximately asymetrical
% hyperboloid characteristic from [t_ini] to [t_fin], using 2*[N]+1
% frequency bands. Extending the time limits will probably work 
% too, as long as [t_ini] is negative and [t_fin] positive. 
% 
% [T1]: backwards temporal masking semiaxis (default= 5)
% [T2]: forwards temporal masking semiaxis (deafult =50)
% [F1]: lower frequencies simultaneous masking semiaxis (default=3)
% [F2]: higher frequencies simulatenous masking semiaxis (default=6)
%
% Author: FVA, 02/2014
error(nargchk(4, 8, nargin));

%% Process input parameters
if (nargin < 5), T1 = 5; end
if (nargin < 6), T2 = 50; end
if (nargin < 7), F1 = 3; end
if (nargin < 8), F2 = 6; end

%% Create the two temporal zones and the two frequency zones
DeltaTBack = t_ini:time_step:0;
DeltaTFor = 0:time_step:t_fin;
DeltaT = [DeltaTBack, DeltaTFor(2:end)];

% The frequency axis is centered around the central analysis frequency.
centralBand = 0;
DeltaFHi = centralBand:N;
DeltaFLo = centralBand-N:centralBand;
DeltaF = [DeltaFLo, DeltaFHi(2:end)];

%Define the hyperboloid: Lower sheet is selected by the minus sign
% a is the semiaxis of left to right: Delta of T (-t__ini, t_fin)
% b is the semiaxis of front to back: Delta of F (positive 0 to N)
% c is the semiaxis of the vertical direction: amount of maskin
hypp = @(x,y,a,b,c)  - (c .* sqrt(1 + x.^2./a^2 + y.^2./b^2));
c = 5;

%% Define four regions
[x,y] = meshgrid(DeltaTBack,DeltaFLo);
SEBackLo = hypp(x,y,T1,F1,c);
[x,y] = meshgrid(DeltaTBack,DeltaFHi);
SEBackHi = hypp(x,y,T1,F2,c);

[x,y] = meshgrid(DeltaTFor,DeltaFLo);
SEForLo = hypp(x,y,T2,F1,c);
[x,y] = meshgrid(DeltaTFor,DeltaFHi);
SEForHi = hypp(x,y,T2,F2,c);

THQ = zeros(length(DeltaF),length(DeltaT));
SE = max(THQ,4*c + ...
    [SEBackLo(2:end,:),SEForLo(2:end,2:end);SEBackHi, SEForHi(:,2:end)]);
% %% Forward amount of masking
% [x,y] = meshgrid(DeltaTFor,DeltaF);
% THQFor = zeros(length(DeltaF),length(DeltaTFor));%Nice way to do it, mainly
% %readjust the forward paraboloid
% 
% SEFor = max(THQFor,4*c + hypp(x,y,50,3,c));
% 
% %% Backward amount of masking
% [x,y] = meshgrid(DeltaTBack,DeltaF);
% THQBack = zeros(length(DeltaF),length(DeltaTBack));%Nice way to do it, mainly
% SEBack = max(THQBack,4*c + hypp(x,y,50,3,c));
% 
%Normalize
SE = (SE - min(min(SE))) ./ (max(max(SE)) - min(min(SE)));
% %PAdding
% addZeros = ((t_fin + t_ini) / time_step) - 1;
% SE = [zeros(N*2 + 1,addZeros),  SEBack, SEFor];
addZeros = ((t_fin + t_ini) / time_step) -1;
SE = [zeros(N*2 +1,addZeros),  SE];
end
