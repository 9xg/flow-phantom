function [outputArg1] = calcSeveringhaus(inputArg1)
%CALCSEVERINGHAUS Summary of this function goes here
%   Detailed explanation goes here
severinghausSO2 = (23400 * (inputArg1.^3 + 150 * inputArg1).^(-1) + 1).^(-1);

outputArg1 = severinghausSO2;
end

