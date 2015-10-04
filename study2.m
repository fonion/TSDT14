%% study2
N = 2^10;
x = randn(N, 1);
n = linspace(0, 1, N);
R0 = 1;
theta_norm = linspace(0,1,N);
theta0 = 0.1; %Dåligt namn ... ? Cutoff?

%% skapar filtrerat brus
[b a] = butter(10, theta0);
filter_noise = filter(b, a, x);


%% Squarer

Ysquare = filter_noise.^2;

%% Half-wave rectifier

Yhalf = stepfunction(filter_noise);

%% AM-SC modulator

omega0 = (2*pi*(theta0/2) + 0.0);

Yamsc = (filter_noise' .* cos(omega0 * n))';

%% Create periodograms

Ysquare_per = fftshift(abs(PeriodFourier(Ysquare)));
Yhalf_per = fftshift(abs(PeriodFourier(Yhalf)));
Yamsc_per = fftshift(abs(PeriodFourier(Yamsc)));

%% Plot estimated PSD:s

figure(1);
subplot(131);
plot(theta_norm,Ysquare_per);
title('Ysquare nonlin')
subplot(132);
plot(theta_norm, Yhalf_per);
title('Yhalf nonlin')
subplot(133);
plot(theta_norm, Yamsc_per);
title('Yamsc nonlin')

%% Histograms
figure(2);
subplot(141)
histogram(Ysquare,100);
title('PSDsquare nonlin')
subplot(142)
histogram(Yhalf,100);
title('PSDhalf nonlin')
subplot(143)
histogram(Yamsc,100);
title('PSDamsc nonlin')
subplot(144)
histogram(filter_noise,100);
title('blackPSD2 lin')

% Varför blir den sista PSD:en inte så gaussisk som vi vill...?
% Kanske så att den är 52 bred och 95 hög...?
% 

%% Nonlinears

figure(3)
plot(thetablack, blackPSD2)
%%
figure(2);
subplot(131)
histoPSDsquare
subplot(132)
histoPSDhalf
subplot(133)
histoPSDamsc

%% Create theoretical PSD:s
%TODO: SKA VI SPEGLA SÅ VI FÅR EN DUBBELSIDIG?
% AM SSB SC? AM SC?
% SKA VI TA BORT DIRAC:EN?
% VARFÖR BLIR VÅR AMSC FÖRSKJUTEN? SKA VI LÄGGA PÅ EN FASFÖRSKJUTNING? 0.2?
%HUR MKT VIKT SKA LÄGGAS PÅ AMPLITUDER NÄR VI JÄMFÖR, TYP ESTIMERADE I LAB1
% AMSC ÄR VÄLDIGT LIK BLACKPSD2 (LINJÄR).
% OM MAN PLOTTAR DEN TEORETISKA RY2 FÅR VI TVÅ STOLPAR. ÄR DET GAUSSISKT?
% KAN EN PSD VARA GAUSSISK?
% VF FÖRSVINNER HALVA PLOTTEN? SAMMA MED VÅRT IDEALA LP.
    % HAR MATLAB PROBLEM MED RECT?

ysquaredC = R0^2/(theta0^2);
ysquared1 = dirac(theta_norm);
ysquared2 = 2/(theta0) * tripuls(theta_norm /theta0);
ysquared3 = 2/(theta0) * tripuls((theta_norm -1 )/theta0);

Ysquared_theor= ysquaredC * (ysquared1 + ysquared2 + ysquared3);

yhalfC = R0/(2 * theta0^2);
yhalf1 = 1;% dirac(theta_norm)/pi;
yhalf2 = 1/2 * rectpuls(theta_norm /theta0);
yhalf3 = theta0/(2 * pi) * tripuls(theta_norm/theta0);
yhalf4 = 1/2 * rectpuls(-theta_norm /theta0);
yhalf5 = theta0/(2 * pi) * tripuls(-theta_norm/theta0);

Yhalf_theor = yhalfC * (yhalf1 + yhalf2 + yhalf3 + yhalf4 + yhalf5);

yamscC = R0/(4 * theta0^2);
yamsc1 = rectpuls((theta_norm+omega0)/theta0);
yamsc2 = rectpuls((theta_norm-omega0)/theta0);
yamsc3 = rectpuls((theta_norm-1+omega0)/theta0);
yamsc4 = rectpuls((theta_norm-1 -omega0)/theta0);

Yamsc_theor = yamscC * (yamsc1 + yamsc2 + yamsc3 + yamsc4);

%% Plot periodograms
figure(2);
subplot(131);
plot(theta_norm,fftshift(abs(Ysquared_theor)));
title('Ysquare nonlin theor')
subplot(132);
plot(theta_norm, fftshift(abs(Yhalf_theor)));
title('Yhalf nonlin theor')
subplot(133);
plot(theta_norm, Yamsc_theor);
title('Yamsc nonlin theor')
%%
figure(3);
Rzhwt = 1/(4*pi)*(tripuls(theta_norm/(2*theta0))+tripuls((theta_norm-1)/(2*theta0)))+...
    +1/4*(rectpuls(theta_norm/(2*theta0))+rectpuls((theta_norm-1)/(2*theta0)));
Rzhwt(1) = Rzhwt(1)+theta0/pi;
plot(theta_norm, fftshift(abs(Rzhwt)));
title('Yamsc nonlin theor david')