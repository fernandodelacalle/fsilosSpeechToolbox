function y=hyperbola(a,b,x)

%Produces the y coordinate of a hyperbola given the x coordinates
%The hyperbola is described by its
% a semimajor axis
% b semiminor axis

y=b-(b*sqrt(1+(x/a).^2));
y=y+b;

end
