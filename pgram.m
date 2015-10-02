function [ periodogram ] = pgram( x )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
N = max(size(x));
Xf = zeros(1,N);

for k = -N/2:N/2-1
    s = 0;
    for n = 0:N-1
        s = s + x(n+1)*exp(-1i*2*pi*k*n/N);
    end
    Xf(k+1+N/2) = s; %The approximate fourier
end
periodogram = 1/(N)*abs(Xf).^2;
end