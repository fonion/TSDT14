%% study2
N = 2^16;
x = randn(N, 1);
n = linspace(0, 1, N);
R0 = 1;
theta_norm = linspace(0,1,N);
theta0 = .1; %D?ligt namn ... ? Cutoff?
fca = theta0/(2*pi);

%% skapar filtrerat brus
[b a] = butter(10, theta0);
filter_noise = filter(b, a, x); %2600

%% Squarer

Ysquare = filter_noise.^2;

%% Half-wave rectifier

Yhalf = stepfunction(filter_noise);

%% AM-SC modulator

omega0 = (2 * pi * fca);
Yamsc = (filter_noise' .* cos(omega0 * n))';

%% Create periodograms

Ysquare_per = abs(pgram(Ysquare));
%Yhalf_per = abs(pgram(Yhalf));
%Yamsc_per = abs(pgram(Yamsc));

%% Plot estimated PSD:s

figure(1);
subplot(131);
plot(theta_norm,Ysquare_per);
%axis([0 1 0 1]);
title('Ysquare')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(132);
plot(theta_norm, Yhalf_per);
%axis([0 1 0 1]);
title('Yhalf')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(133);
plot(theta_norm, Yamsc_per);
%axis([0 1 0 1]);
title('Yamsc')
xlabel('Theta')
ylabel('Power Spectral Density')

%% Histograms
figure(2);

subplot(141)
histogram(Ysquare,100);
axis([-2 2 0 25000])
title('square')
xlabel('amplitude')
ylabel('number of samples')

subplot(142)
histogram(Yhalf,100);
axis([-2 2 0 35000])
title('half')
xlabel('amplitude')
ylabel('number of samples')



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
% 

%% Nonlinears

%figure(3)
%plot(thetablack, blackPSD2)
%%
%figure(2);
%subplot(131)
%histoPSDsquare
%subplot(132)
%histoPSDhalf
%subplot(133)
%histoPSDamsc

%% Create theoretical PSD:s
ysquaredC = R0^2/(theta0^2);
ysquared1 = 0;%dirac(theta_norm);
ysquared2 = 2/(theta0) * tripuls(theta_norm /theta0);
ysquared3 = 2/(theta0) * tripuls((theta_norm -1 )/theta0);

Ysquared_theor= ysquaredC * (ysquared1 + ysquared2 + ysquared3);

yhalfC = R0/(2 * theta0^2);
yhalf1 = 0;% dirac(theta_norm)/pi;
yhalf2 = 1/2 * rectpuls(theta_norm /theta0);
yhalf3 = theta0/(1 * pi) * tripuls(theta_norm/theta0);
yhalf4 = 1/2 * rectpuls((theta_norm -1)/theta0);
yhalf5 = theta0/(1 * pi) * tripuls((theta_norm - 1 )/theta0);

Yhalf_theor = yhalfC * (yhalf1 + yhalf2 + yhalf3 + yhalf4 + yhalf5);

yamscC = R0/(4 * theta0^2);
yamsc1 = rectpuls((theta_norm + theta0 * omega0)/theta0);
yamsc2 = rectpuls((theta_norm - theta0 * omega0)/theta0);
yamsc3 = rectpuls((theta_norm - 1 + theta0 * omega0)/theta0);
yamsc4 = rectpuls((theta_norm - 1 - theta0 * omega0)/theta0);

Yamsc_theor = yamscC * (yamsc1 + yamsc2 + yamsc3 + yamsc4);

%% Plot periodograms
figure(2);
subplot(131);
plot(theta_norm,Ysquared_theor);
title('Ysquare nonlin theor')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(132);
plot(theta_norm, Yhalf_theor);
title('Yhalf nonlin theor')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(133);
plot(theta_norm, Yamsc_theor);
title('Yamsc nonlin theor')
xlabel('Theta')
ylabel('Power Spectral Density')

%% asdf
% NÄR VI RAPPORTERAR: SÄG ATT DIRAC:EN I 0:AN FÅR VÅR BILD ATT SE SKEV UT.
% PLOTTA T.EX. FRÅN 0->1 I HÖJD. DÅ ÄR DEN JÄTTE LIK IDEALA xDxDxDxDxD


%%
%figure(3);
%Rzhwt = 1/(4*pi)*(tripuls(theta_norm/(2*theta0))+tripuls((theta_norm-1)/(2*theta0)))+...
%    +1/4*(rectpuls(theta_norm/(2*theta0))+rectpuls((theta_norm-1)/(2*theta0)));
%Rzhwt(1) = Rzhwt(1)+theta0/pi;
%plot(theta_norm, fftshift(abs(Rzhwt)));
%title('Yamsc nonlin theor david')


%% NYA FRÅGOR
% 1. VI TROR ATT ALLA THETA0 I DETTA SCRIPT SKA VA SAMMA. STÄMMER DET?



%% FR?GOR
% 1. VARF?R BLIR V?R AMSC F?RSKJUTEN? 
%       SKA VI L?GGA P? EN FASF?RSKJUTNING? 0.2?
% 2. VI BORDE HA H?G AMPLITUD I NOLL. O?NDLIGT 
%       ANTAL SAMPELS --> O?NDLIG H?JD I NOLL. 
%       TY DET ?R D?R DIRAC:EN H?RJAR
% 3. HISTO: AMSC ?R V?LDIGT LIK BLACKPSD2 (LINJ?R).
        %VAD ?R SKILLNADEN?
% 4. KAN EN PSD VARA GAUSSISK?
% 5. Why do we plot the periodograms? 
% 6. In order to see if there are any new frequencies.
%       If so - the system is not LTI.
% 7. Why do we plot the histograms? 
%       In order to see that out signal is gaussian. 
%       If not - the system is not LTI.
% 6. ASDF
