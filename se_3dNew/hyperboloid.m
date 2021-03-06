%% drawing a right hyperboloid of two sheets as per:
% http://www.ehow.com/how_8716313_plot-hyperboloid-matlab.html

% 1.- Solve for z in the hyperboloid of two sheets
% 2.- write a Matlab function for it
a = 2;%x semixis
b = 3;%y semiaxis
c = 5;%z semiaxis
hyp = @(x,y)  - (c .* sqrt(1 + x.^2./a^2 + y.^2./b^2));

%% 3 DRAW with ezsurf
ezsurf(hyp, [-10, 20, -15, 20])


%% 4 Now parametrize in the semi-axes
hypp = @(x,y,a,b,c)  - (c .* sqrt(1 + x.^2./a^2 + y.^2./b^2));



%% 5 Draw several hyperboloids
ezsurf(@(x,y) hypp(x, y, 2, 3, 4), [-1 10 -20 20])

hold on

ezsurf(@(x,y) hypp(x, y, 1,3, 4), [-1 10 -20 20])


% NOw rather use surf
[x,y] = meshgrid(linspace(-1,10,60),linspace(-20,20,60));
z = hypp(x,y,1,3,4);
surf(x,y,z)
xlabel('Delta T')
ylabel('Delta F')
zlabel('Masking')

%OK So now go back to the end of creaSE3dNew to see this in
% (DeltaF,DeltaT,M) coordinates
