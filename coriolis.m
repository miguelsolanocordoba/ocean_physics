function f = coriolis(lat)
%f = coriolis(lat) computes Coriolis frequency. 
% 
% Created: Solano, M. July 7, 2020

sig = 7.292115E-5;
f = 2*sig.*sin(lat); 
