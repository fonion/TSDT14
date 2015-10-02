%% generate noise

 

x = randn(2^10, 1);

 

%% filter low degree

 

R0 = 1;
a = 0.9;
N = 2^10;


theta = linspace(0,1,N);

H1 = 1./(1-a*exp(-1i*2*pi.*theta)); 


%PSD low degree

Ry1 = fftshift(R0./(1+a^2-2*a.*cos(2*pi.*theta)));

plot(theta, Ry1);
title('PSD low degree filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%%

 

%ACF low degree

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

%plot(theta, H2)

 

%PSD idealt

Ry2 = R0/(theta0^2).*rectangularPulse(theta/theta0);

plot(theta, Ry2)
title('PSD idealt filter')
xlabel('Theta')
ylabel('Power Spectral Density')
 %%

%ACF idealt

n1 = linspace((-N/2), (N/2), N+1);
n22 = linspace(-20, 20, 100);

ry21 = R0/(theta0)*sinc(theta0*n1); %intervall alla n
ry22 = R0/(theta0)*sinc(theta0*n22); %intervall -20 till 20

 

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

 

 

 

%% Estimation of PSD and ACF using Bartlett
% skapar h1, h2 och saker att plotta mot
 

h1 = ifft(H1);
h2 = ifft(H2);


n3 = linspace((-N), (N), (2*N)+1);
n4 = linspace((-20), 20, (N));

theta2 = linspace(0,1,(2*N)+1);
theta3 = linspace(0, 1, 2058); 
 
%% outputs

%outputh1 = real(conv(x, h1));

bf = 1;
af = [1 ;-a];
outputh1 = filter(bf,af,x);
outputh2 = real(conv(x, h2));

%% test nya
ACFH1 = ACF_estimation(outputh1, 'bartlett');

%% Estimations

acfH1 = (ACF_bart(outputh1, N));
psdH1 = fftshift(abs(fft(acfH1)));

acfH2 = (ACF_bart(outputh2, N));
psdH2 = fftshift(abs(fft(acfH2)));

%% Blackman, smoothing osv

W = blackman(10);

FW = real(fft(W));

%WpsdH1 = abs(conv(FW, psdH1)); 
%WpsdH2 = abs(conv(FW, psdH2));

WacfH1= abs(conv(FW, acfH1));

blackPSDH1 = (fftshift(abs(fft(WacfH1)))) ;

%% plots

figure(2)
subplot(221)
plot(n3, acfH1) 
title ('ACF H1');

subplot(222)
plot(theta3, blackPSDH1)
title ('PSD H1 med blackman');

subplot(223)
plot(n3, acfH2)
title ('ACF H2');

subplot(224)
plot(theta3, WpsdH2)
title ('PSD H2 med blackman');

%% plottar mot varandra. 
% framräknade ACF och estimerade ACFH1

subplot(121)
stem(n2, ry12)
title('uträknade ACF för H1')


subplot(122)
stem(n4, ACFH1) 
title ('estimerade ACF för H1');

%% framräknade ACF och estimerade ACFH2

figure(2)
subplot(121)
stem(n22, ry22)
title('uträknade ACF H2')


subplot(122)
stem(n4, acfH2) 
title ('estimerade ACF för H2');



