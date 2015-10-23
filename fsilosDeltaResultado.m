function [ delta ] = deltaResultado( p )
%DELTARESULTADO Summary of this function goes here
%   Detailed explanation goes here
n = 1560*4;
delta = 1.95* sqrt((p*(100-p))/(n));

end

