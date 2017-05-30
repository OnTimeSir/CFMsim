function [ vout ] = Velocity8( v )
%VELOCITY8 Summary of this function goes here
%   Detailed explanation goes here
vout = int8(mod(double(v)+128,256)-128);
end

