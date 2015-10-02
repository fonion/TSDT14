function [ACF_hat] = ACF_help(y,i)
    N=length(y);
    ACF_hat = 0;
    for n = 1:N-abs(i)-1
        ACF_hat = ACF_hat+(y(n+abs(i))*y(n))/N;
    end
end