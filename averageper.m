function [ PSDaver ] = averageper(x, k)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
N = max(size(x));
Nk = N/k;
%rhat = zeros(1,N);
rhat = 0;
%size(rhat)
%size(ACF_estimation(x(i:Nk),'Bartlett'))
for i = 1:Nk:N
    rhat = [rhat , ACF_estimation(x(i:Nk),'Bartlett')]; %get elements
             %check ACF
             %sum all acf, divide by amount
end

rhat = rhat/k;
PSDaver = fft(rhat);

end

