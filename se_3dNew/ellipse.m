function y=ellipse(a,b,x)

%Produces the y coordinate of an ellipse given the x coordinates
%The ellipse is described by its
% a semimajor axis
% b semiminor axis

y=(b*sqrt(1-(x/a).^2));

end
