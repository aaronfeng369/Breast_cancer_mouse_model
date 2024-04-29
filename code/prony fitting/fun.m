function [ MSE ] = fun( x,t,data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    G=x(1)-(x(2)*(1-exp((-1*t)/x(5)))+x(3)*(1-exp((-1*t)/x(6)))+x(4)*(1-exp((-1*t)/x(7))));
    MSE=sum((data-G).^2)/14000;


end
