%% study2
N = 2^10;
x = randn(N, 1);
n = linspace(0, N, N);

%% skapar filtrerat brus

[b a] = butter(10, 0.2);
filter_noise = filter(b, a, x);


%% Squarer

Ysquare = filter_noise.^2;

%% Half-wave rectifier

Yhalf = stepfunction(filter_noise);

%% AM-SC modulator
cutoff = 0.1;
theta0 = (2*pi*(cutoff));

Yamsc = filter_noise* cos(theta0 * n);

%% 

