function [outputArg1] = calcKelman(inputArg1,inputArg2)
% po2, temp
%CALCKELMAN Summary of this function goes here
%   Detailed explanation goes here
% 1st is pO2, 2nd is temp

a1 = -8.5322289e3; 
a2 = 2.1214010e3;
a3 = -6.7073989e1;
a4 = 9.3596087e5;
a5 = -3.1346258e4;
a6 = 2.3961674e3;
a7 = -6.7104406e1;

xKelman = inputArg1.*10.^(0.0024*(37-inputArg2)+0.4*(7.40-7.40)+ 0.06*(log(40)-log(40)));
kelmanSO2 = ((a1.*xKelman+a2*xKelman.^2+a3*xKelman.^3+xKelman.^4) ./ (a4+a5*xKelman+a6*xKelman.^2 +a7*xKelman.^3 + xKelman.^4));

outputArg1 = kelmanSO2;
end

