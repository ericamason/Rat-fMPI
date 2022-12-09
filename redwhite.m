function c = redwhite
%REDWHITE    Shades of red and white color map

r = 256*ones(256,1)/256;
g = [256:-1:1]'/256;
b = g;

c = [r g b]; 
