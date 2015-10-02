function [ ry_vector ] = ACF_estimation(x,type)
N = length(x);
ry_vector = zeros(1,N);

for k = -N/2+1:N/2
    sum = 0;
    for n = 1:(N-abs(k))
        sum = sum+x(n+abs(k))*x(n);
    end
    ry_vector(k+N/2) = sum;
    switch type
    case 'Blackman'
        ry_vector(k+N/2) = ry_vector(k+N/2)/(N-abs(k));
    case 'Bartlett'
        ry_vector(k+N/2) = ry_vector(k+N/2)/(N); 
    end
end     
end