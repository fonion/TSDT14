%% generate noise

 

x = randn(2^10, 1);

 

%% filter low degree

 

R0 = 1;

a = 0.9;

N = 2^10;

 

 

 

theta = linspace(0,1,N);

 

H1 = 1./(1-a*exp(-1i*2*pi.*theta)); 

 

%PSD

Ry1 = R0./(1+a^2-2*a.*cos(2*pi.*theta));

 

plot(theta, Ry1);

%%

 

%ACF

n1 = linspace((-N/2), (N/2), N+1);

n2 = linspace(-20, 20, 41);

 

ry11 = (R0/(1-a^2)).*a.^abs(n1); %intervall alla n

ry12 = (R0/(1-a^2)).*a.^abs(n2); %intervall -20 till 20

 

subplot(121)

plot(n1, ry11)

subplot(122)

stem(n2, ry12)

 

%% Idealt filter

 

theta0 = 0.7;

 

H2 = (1/theta0)*rectangularPulse(theta/theta0);

%plot(theta, H2)

 

%PSD

Ry2 = R0/(theta0^2).*rectangularPulse(theta/theta0);

plot(theta, Ry2)

 

%ACF

n1 = linspace((-N/2), (N/2), N+1);

n2 = linspace(-20, 20, 100);

 

ry21 = R0/(theta0)*sinc(theta0*n1); %intervall alla n

ry22 = R0/(theta0)*sinc(theta0*n2); %intervall -20 till 20

 

subplot(121)

plot(n1, ry21)

subplot(122)

stem(n2, ry22)

 

 

 

%% Estimation of PSD and ACF using Bartlett

 

h1 = ifft(H1);
h2 = ifft(H2);


n3 = linspace((-N), (N), (2*N)+1);

theta2 = linspace(0,1,(2*N)+1);

 

%outputs

outputh1 = real(conv(x, h1));
outputh2 = real(conv(x, h2));

acfH1 = abs(ACF_bart(outputh1, N));
psdH1 = fftshift(abs(fft(acfH1)));

acfH2 = abs(ACF_bart(outputh2, N));
psdH2 = fftshift(abs(fft(acfH2)));

 

 

subplot(231)
title('ACF H1');
plot(n3, acfH1) 

subplot(232)
plot(theta2, psdH1)

subplot(233)
plot(n3, acfH2)

subplot(234)
plot(theta2, psdH2)

