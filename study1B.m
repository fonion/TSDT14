%% generate noise

x = randn(2^16, 1);
figure(100)
%% filter low degree

R0 = 1;
a = 0.9;
N = 2^16;

theta = linspace(0,1,N);

H1 = 1 ./ (1 - a* exp(-1i * 2 * pi .* theta));

%% PSD low degree

Ry1 = (R0./(1+a^2-2*a.*cos(2*pi.*theta)));
%%
figure(8)
plot(theta, Ry1);
title('PSD low degree filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% ACF low degree

n1 = linspace((-N/2), (N/2), N+1);
n2 = linspace(-20, 20, 41);

ry11 = (R0/(1-a^2)).*a.^abs(n1); %intervall alla n
ry12 = (R0/(1-a^2)).*a.^abs(n2); %intervall -20 till 20

figure(2)
subplot(121)
plot(n1, ry11)
title('ACF lowdegree filter f?r alla n')
xlabel('sampels')
ylabel('Auto Correlation Function')

subplot(122)
stem(n2, ry12)
title('ACF lowdegree filter -20 <n< 20')
xlabel('sampels')
ylabel('Auto Correlation Function')

 

%% Idealt filter

thetaideal = linspace(0, 1, N); 

theta0 = 0.1;
H2 = rectangularPulse(thetaideal/(2*theta0));

%plot(theta, H2)

%% PSD idealt

Ry2 = R0*rectangularPulse(thetaideal/(2*theta0));

Ry2 = Ry2 + R0*rectangularPulse((thetaideal - 1)/(2*theta0));
%%
figure(2)
plot(thetaideal, Ry2)
axis([0 1 0 1.5])
title('PSD idealt filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% ACF idealt

n1 = linspace((-N/2), (N/2), N+1);
n22 = linspace(-20, 20, 41);

ry21 = R0*theta0*2*sinc(2*theta0*n1); %intervall alla n
ry22 = R0*theta0*2*sinc(2*theta0*n22); %intervall -20 till 20

 
figure(1)
subplot(121)
plot(n1, ry21)
title ('ACF Idealt filter alla n');
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

ACFh1 = ACF_estimation(outputh1, 'Bartlett'); %OBS case-sensitive...

%% Plot ACF Low degree raw

twenty = (length(ACFh1)/2)-20;

figure(1)
subplot(121)
plot(n3, ACFh1)
title('raw estimate ACF, all n')
xlabel('samples')
ylabel('Auto Correlation Function')

subplot(122)
stem(n4, ACFh1(twenty:(twenty+40)))
title('raw estimate ACF, -20 < n < 20')
xlabel('samples')
ylabel('Auto Correlation Function')

%% Estimate PSD low degree raw

PSDh1 = abs((fft(ACFh1)));

figure(3)
plot(theta, PSDh1);
title('raw estimate PSD low degree filter')
xlabel('Theta')
ylabel('Power Spectral Density')
65537
11
%%

window = blackman(11);
N = 65537-11
zeropaddedwindow = [zeros(1,(N)/2) window zeros(1,(N)/2)];
length(zeropaddedwindow)
%%

%% Periodogram = PSD f?r low degree.
figure(4)
PSDgram = fftshift(pgram(ACFh1));

plot(theta, PSDgram)
title('periodogram estimate PSD low degree filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% averaging på periodogram PSD aver1
intervals = 2^8;

PSDaver1 = averageper(outputh1,intervals);

figure(3);
thetaaver = linspace(0,1,N/intervals + 1);
plot(thetaaver, (abs(PSDaver1)));
axis([0 1 0 1]);
title('PSD averaging low degree filter')
xlabel('Theta')
ylabel('Power Spectral Density')
%% smoothing p? averaging PSD averwindow!

ACFaver1 = ifft(PSDaver1);
window = blackman(length(ACFaver1))';
ACFaver1window = window .* ACFaver1;
PSDaver1window = abs(fft (ACFaver1window));

figure(4)
plot(thetaaver, PSDaver1window);
axis([0 1 0 1]);
title('PSD averaging and smoothing low degree filter')
xlabel('Theta')
ylabel('Power Spectral Density')
%% averaging ACF1
averinterval = 1;
avertemp = averageper(outputh1,averinterval);
ACFaver1 = ifft(avertemp);
n5 = linspace(-N/2 , N/2 , max(size(ACFaver1)));

twenty1 = (length(ACFaver1)/2)-19;
figure(6)
subplot(121)
plot(n5, ACFaver1)
title('averaging ACF, all n')
xlabel('sampels')
ylabel('Auto Correlation Function')

subplot(122)
stem(n4, ACFaver1(twenty1:(twenty1+40)))
title('averaging ACF, -20 < n < 20')
xlabel('sampels')
ylabel('Auto Correlation Function')

%% smoothing p? aver ACF1


%window = blackman(length(ACFaver1'));
windowlength =  11;
window = blackman(windowlength )';
padding = length(ACFaver1)-windowlength ;
window = [zeros(1,padding/2) window zeros(1,padding/2)];
ACFaver1window = window .* ACFaver1;

n6 = linspace(-N/2, N/2, max(size(window))); 
twenty2 = (length(ACFaver1window)/2)-20;

figure(3)
subplot(121)
plot(n6, ACFaver1window)
title('Smoothing, all n')
xlabel('samples')
ylabel('Auto Correlation Function')

subplot(122)
stem(n4, ACFaver1(twenty2:(twenty2+40)))
title('Smoothing, -20 < n < 20')
xlabel('samples')
ylabel('Auto Correlation Function')

%% 
%
%
% Ideal filter
%
%

%%


[b a] = butter(10, 2*theta0);
outputh2 = filter(b, a, x);


ACFh2 = ACF_estimation(outputh2, 'Bartlett');

%% Plot Ideal

figure(3)
subplot(121)
plot(n3, ACFh2)
title('raw estimate ACF, all n')
xlabel('samples')
ylabel('Auto Correlation Function')


twenty3 = (length(ACFh2)/2)-20;
subplot(122)
stem(n4, ACFh2(twenty3:(twenty3+40)))
title('raw estimate ACF, -20 < n < 20')
xlabel('samples')
ylabel('Auto Correlation Function')


%% raw PSD
PSDraw2 = abs(fft(ACFh2));
plot(theta, PSDraw2)
title('raw PSD, ideal filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% Periodogram = PSD f?r Ideal. raw estimate

PSDgram2 = fftshift(pgram(ACFh2));

figure(1)
plot(theta, PSDgram2)
title('PSD, periodogram ideal filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% Smoothing av estimated PSD Blackman

% windowsize = 3;
% blackmanw = blackman(windowsize);
% %blackACF = BlackmanW*ACFh1;
% BlackmanW = fft(blackmanw);
% blackPSD2 = abs(conv(BlackmanW, PSDgram2));
% thetablack = linspace(0, 1, 1024+windowsize-1);
% 
% figure(2)
% plot(thetablack, blackPSD2);

%% averaging 



intervals = 2^6;

thetaaver = linspace(0,1,N/intervals + 1);
PSDaver2 = averageper(outputh2,intervals);

figure(1);
plot(thetaaver, (abs(PSDaver2)));
title('PSD, averaging ideal filter')
xlabel('Theta')
ylabel('Power Spectral Density')

%% averaging by johan

interval = 2^6;
[plotmatrix, aver2] = averaging( PSDraw2, interval );
thetaaver = linspace(0,1,length(aver2));

figure(1);
plot(plotmatrix, aver2);
title('PSD, averaging ideal filter')
xlabel('Theta')
ylabel('Power Spectral Density')



%% smoothing på averaging
windowlength = 513;
window = blackman(windowlength)/(windowlength/2);
% padding = length(ACFaver2)-windowlength ;
% window = [zeros(1,padding/2) window zeros(1,padding/2)];
% 
% ACFaver2window = window .* ACFaver2;
% PSDaver2window = abs(fft (ACFaver2window));

PSDaver2window = filter(window,1,PSDraw2);

%%fixa att plotta mot rätt theta!! :D

figure(2)
plot(PSDaver2window);
%axis([0 1 0 3]);
%% smoothing p? averaging

ACFaver2 = ifft(PSDaver2);
windowlength = 49;
window = blackman(windowlength)';
padding = length(ACFaver2) - windowlength;
window = [zeros(1,padding/2) window zeros(1,padding/2)];
ACFaver2window = window .* ACFaver2;
PSDaver2window = abs(fft (ACFaver2window));

figure(2)
plot(thetaaver, PSDaver2window);
%axis([0 1 0 50]);
>>>>>>> d6e12f952663d6320b29d45ffa4eb0e8e647036a
title('PSD, averaging and smoothing')
xlabel('Theta')
ylabel('Power Spectral Density')
ylim([0,6])
%% averaging ACF2
averinterval = 1;
avertemp = averageper(outputh2,averinterval);
ACFaver2 = ifft(avertemp);
n5 = linspace(-N/2 , N/2 , max(size(ACFaver2)));

twenty4 = (length(ACFaver2)/2)-20;

figure(7)
subplot(121)
plot(n5, ACFaver2)
title('averaging ACF2, all n')
xlabel('sampels')
ylabel('Auto Correlation Function')

subplot(122)
stem(n4, ACFaver2(twenty4:(twenty4+40)))
title('averaging ACF, -20 < n < 20')
xlabel('sampels')
ylabel('Auto Correlation Function')

% subplot(122)
% stem(n4, ACFaver1)
% title('averaging ACF, -20 < n < 20')


%% smoothing på aver ACF2
%window = blackman(length(ACFaver2'));

windowlength =  101;
window = blackman(windowlength )';
padding = length(ACFaver2)-windowlength ;
window = [zeros(1,padding/2) window zeros(1,padding/2)];

ACFaver2window = window .* ACFaver2;
n6 = linspace(-N/2, N/2, max(size(window)));

twenty5 = (length(ACFaver2window)/2)-20;

figure(7)
subplot(121)
plot(n6, ACFaver2window)
title('smoothing, all n')
xlabel('samples')
ylabel('Auto Correlation Function')

subplot(122)
stem(n4, ACFaver2window(twenty5:(twenty5+40)))
title('ACF, smoothing, -20 < n < 20')
xlabel('samples')
ylabel('Auto Correlation Function')
 