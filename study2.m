%% study2
N = 2^18;
x = randn(N, 1);
n = linspace(0, 1, N);
R0 = 1;
theta_norm = linspace(0,1,N);
theta0 = .05; %D?ligt namn ... ? Cutoff?
fca = theta0*8;
omega0 = (2 * pi * fca);

%% skapar filtrerat brus
[b a] = butter(10, 2 * theta0);
filter_noise = filter(b, a, x); %2600

%% Squarer

Ysquare = filter_noise.^2;

%% Half-wave rectifier

Yhalf = stepfunction(filter_noise);

%% AM-SC modulator

Yamsc = (filter_noise' .* cos(omega0.* theta_norm))';

%% Create periodograms

Ysquare_per = 1/(length(Ysquare))*abs(fft(Ysquare)).^2;%abs(pgram(Ysquare));
Yhalf_per = 1/(length(Yhalf))*abs(fft(Yhalf)).^2;%abs(pgram(Yhalf));
Yamsc_per = 1/(length(Yamsc))*abs(fft(Yamsc)).^2;%abs(pgram(Yamsc));


%% Plot estimated PSD:s

figure(1);
subplot(131);
plot(theta_norm,(Ysquare_per));
axis([0 1 0 4]);
title('Ysquare')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(132);
plot(theta_norm, (Yhalf_per));
axis([0 0.25 0 4]);
title('Yhalf')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(133);
plot(theta_norm, (Yamsc_per));
axis([0 1 0 4]);
title('Yamsc')
xlabel('Theta')
ylabel('Power Spectral Density')


%% Histograms
figure(2);

subplot(141)
histogram(Ysquare,100);
axis([-5 5 0 100000])
title('square')
xlabel('amplitude')
ylabel('number of samples')

subplot(142)
histogram(Yhalf,100);
axis([-5 5 0 150000])
title('half')
xlabel('amplitude')
ylabel('number of samples')

title('Yhalf nonlin theor')


subplot(143)
histogram(Yamsc,100);
title('amsc')
xlabel('amplitude')
ylabel('number of samples')


subplot(144)
histogram(filter_noise,100);
title('gaussian')
xlabel('amplitude')
ylabel('number of samples')


% Varf?r blir den sista PSD:en inte s? gaussisk som vi vill...?
% Kanske s? att den ?r 52 bred och 95 h?g...?

%% Create theoretical PSD:s
ysquaredC = R0^2;% R0^
ysquared1 = 0;%dirac(theta_norm);
ysquared2 = 2 * tripuls(theta_norm /(4*theta0));
ysquared3 = 2 * tripuls((theta_norm-1) /(4*theta0));

Ysquared_theor= ysquaredC * (ysquared1 + ysquared2 + ysquared3);

yhalfC = R0;
yhalf1 = 0;% dirac(theta_norm)/pi;
yhalf2 = 1/4 * rectpuls(theta_norm /(2*theta0));
yhalf3 = 1/(4 * pi * theta0) * tripuls(theta_norm/(4*theta0));
yhalf4 = 1/4 * rectpuls((theta_norm -1)/(2*theta0));
yhalf5 = 1/(4 * pi * theta0) * tripuls((theta_norm -1)/(4*theta0));

Yhalf_theor = yhalfC * (yhalf1 + yhalf2 + yhalf3 + yhalf4 + yhalf5);

yamscC = R0/4;
yamsc1 = rectpuls((theta_norm + theta0 * omega0)/(2*theta0));
yamsc2 = rectpuls((theta_norm - theta0 * omega0)/(2*theta0));
yamsc3 = rectpuls((theta_norm - 1 + theta0 * omega0)/(2*theta0));
yamsc4 = rectpuls((theta_norm - 1 - theta0 * omega0)/(2*theta0));

Yamsc_theor = yamscC * (yamsc1 + yamsc2 + yamsc3 + yamsc4);

%% Plot periodograms
figure(3);
subplot(131);
plot(theta_norm,Ysquared_theor);
axis([0 0.25 0 5]);
title('Ysquare nonlin theor')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(132);
plot(theta_norm, Yhalf_theor);
title('Yhalf nonlin theor')
xlabel('Theta')
ylabel('Power Spectral Density')
axis([0 0.25 0 5]);
subplot(133);
plot(theta_norm, Yamsc_theor);
title('Yamsc nonlin theor')
xlabel('Theta')
ylabel('Power Spectral Density')
axis([0 1 0 2]);