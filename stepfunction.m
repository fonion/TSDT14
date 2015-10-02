function [ y ] = stepfunction( x )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
N = size(x);
y = zeros(N);
for i = 1:N
    
    if x(i) > 0
        y(i) = x(i);
    else
        y(i) = 0;
    end
end

end

