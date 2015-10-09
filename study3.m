% Study 3

%% skapar filtrerat brus
N = 2^10;
x = randn(1, N);

[b a] = butter(10, theta0);
filter_noise = filter(b, a, x);

%% 

% +1 -1

n7 = linspace(0, N, N);
Ypm = filter_noise .* (-1) .^ n7;

% 010101

Y01 = (filter_noise + filter_noise .* (-1) .^ n7)/2;

%% estimate PSD (periodogram)

Ypm_per = abs(PeriodFourier(Ypm));
Y01_per = abs(PeriodFourier(Y01));

figure(1);
subplot(121);
plot(theta_norm,Ypm_per);
axis([0 1 0 1]);
title('Ysquare nonlin')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(122);
plot(theta_norm, Y01_per);
axis([0 1 0 1]);
title('Yhalf nonlin')
xlabel('Theta')
ylabel('Power Spectral Density')

