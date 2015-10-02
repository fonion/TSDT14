%% generate noise

x = randn(2^10, 1);

%% filter low degree

R0 = 1;
a = 0.9;
N = 2^10;

theta = linspace(0,1,N);

H1 = 1 ./ (1 - a* exp(-1i * 2 * pi .* theta));

%% PSD low degree

Ry1 = fftshift(R0./(1+a^2-2*a.*cos(2*pi.*theta)));
%%
figure(1)
plot(theta, Ry1);
title('PSD low degree filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% ACF low degree

n1 = linspace((-N/2), (N/2), N+1);
n2 = linspace(-20, 20, 41);

ry11 = (R0/(1-a^2)).*a.^abs(n1); %intervall alla n
ry12 = (R0/(1-a^2)).*a.^abs(n2); %intervall -20 till 20


subplot(121)
plot(n1, ry11)
title('ACF lowdegree filter för alla n')
xlabel('sampels')
ylabel('Auto Correlation Function')

subplot(122)
stem(n2, ry12)
title('ACF lowdegree filter -20 <n< 20')
xlabel('sampels')
ylabel('Auto Correlation Function')

 

%% Idealt filter

 

theta0 = 0.1;
H2 = (1/theta0)*rectangularPulse(theta/theta0);

plot(theta, H2)

%% PSD idealt

Ry2 = R0/(theta0^2).*rectangularPulse(theta/theta0);

Ry2 = fftshift(Ry2);

figure(2)
plot(theta, Ry2)
title('PSD idealt filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% ACF idealt

n1 = linspace((-N/2), (N/2), N+1);
n22 = linspace(-20, 20, 41);

ry21 = R0/(theta0)*sinc(theta0*n1); %intervall alla n
ry22 = R0/(theta0)*sinc(theta0*n22); %intervall -20 till 20

 
figure(1)
subplot(121)
plot(n1, ry21)
title ('ACF Idealt filter alla n');
title('ACF lowdegree filter för alla n')
xlabel('sampels')
ylabel('Auto Correlation Function')


subplot(122)
stem(n22, ry22)
title ('ACF Idealt filter -20 < n < 20');
xlabel('sampels')
ylabel('Auto Correlation Function')

%% ESTIMATION
%
%
%

%% 
%
%Low degree fiter 
%
%
%%
n3 = linspace((-N/2), (N/2), N);
n4 = linspace(-20, 20, 41);
bf = 1;
af = [1;-a];
outputh1 = filter(bf,af,x);

ACFh1 = ACF_estimation(outputh1, 'bartlett');
%% Plot ACF Low degree
figure(2)
subplot(121)
plot(n3, ACFh1)
title('estimate ACF, all n')


subplot(122)
stem(n4, ACFh1(492:(492+40)))
title('estimate ACF, -20 < n < 20')

%% Estimate PSD low degree
figure(2)
PSDh1 = abs(fftshift(fft(ACFh1)));

figure(3)
plot(theta, PSDh1);
title('PSD low degree filter')

%% Smoothing av estimated PSD Blackman

windowsize = 3;
blackmanw = blackman(windowsize);
%blackACF = BlackmanW*ACFh1;
BlackmanW = fft(blackmanw);
blackPSD = abs(conv(BlackmanW, PSDh1));
thetablack = linspace(0, 1, 1024+windowsize-1);

figure(2)
plot(thetablack, blackPSD);

%% Periodogram = PSD för low degree. raw estimate

PSDgram = pgram(ACFh1);

plot(theta, PSDgram)

%% averaging på periodogram
intervals = 45;
PSDaver1 = averageper(outputh1,intervals);
figure(3);
thetaaver = linspace(0,1,N/intervals + 1);
plot(thetaaver, fftshift(abs(PSDaver1)));
title('PSDaver1')


%% 
%
%
% Ideal filter
%
%

%%

[b a] = butter(10, 0.1);
outputh2 = filter(b, a, x);

ACFh2 = ACF_estimation(outputh2, 'bartlett');

%% Plot Ideal

figure(2)
subplot(121)
plot(n3, ACFh2)
title('estimate ACF, all n')


subplot(122)
stem(n4, ACFh2(492:(492+40)))
title('estimate ACF, -20 < n < 20')

%% Periodogram = PSD för Ideal. raw estimate

PSDgram2 = pgram(ACFh2);

figure(1)
plot(theta, PSDgram2)


%% Smoothing av estimated PSD Blackman

windowsize = 3;
blackmanw = blackman(windowsize);
%blackACF = BlackmanW*ACFh1;
BlackmanW = fft(blackmanw);
blackPSD2 = abs(conv(BlackmanW, PSDgram2));
thetablack = linspace(0, 1, 1024+windowsize-1);

figure(2)
plot(thetablack, blackPSD2);

%% averaging på periodogram
thetaaver = linspace(0,1,N/intervals + 1);

intervals = 45;
PSDaver2 = averageper(outputh2,intervals);

figure(1);
plot(thetaaver, fftshift(abs(PSDaver2)));
title('PSDaver2')
 