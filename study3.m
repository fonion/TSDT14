% Study 3

N = 2^10;
x = randn(N, 1);
R0 = 1;
theta_norm = linspace(0,1,N);
theta0 = .1; %D?ligt namn ... ? Cutoff?
fca = theta0/(2*pi);

%% skapar filtrerat brus
>>>>>>> d6e12f952663d6320b29d45ffa4eb0e8e647036a

[b a] = butter(10, theta0);
filter_noise = filter(4.6*b, a, x);

%%
% onematrix = ones(1, N);
% 
% for i:N
%     if mod(i,2)
%         i = 0;
%     end 
% end
%% 

% +1 -1

n7 = linspace(0, N, N);
Ypm = filter_noise .* (-1) .^ n7;

% 010101

phi = 1; %randi([0 1], 1, N);

Y01 = (filter_noise - filter_noise .* (-1) .^ (n7 .* phi))/2;

%% estimate PSD (periodogram)

Ypm_per = (abs(pgram(Ypm)));
Y01_per = (abs(pgram(Y01)));

%%
Ypm_per = ifftshift(Ypm_per);


figure(1);
subplot(121);
plot(theta_norm,Ypm_per);
%axis([0 1 0 40]);
title('PSD, Ypm')
xlabel('Theta')
ylabel('Power Spectral Density')
subplot(122);
plot(theta_norm, Y01_per);


axis([0 1 0 10]);

title('PSD, Y01')
xlabel('Theta')
ylabel('Power Spectral Density')

%% Create theoretical PSD:s

Rypm = (R0/(theta0)^2) .* rectpuls((theta_norm-0.5)/theta0);

figure(3)
plot(theta_norm, Rypm)
title('Theoretical PSD')
xlabel('theta')
ylabel('Power Spectral Density')

Ry011 = (R0/(4*theta0^2)) .* rectpuls(theta_norm/theta0);
Ry0111 = (R0/(4*theta0^2)) .* rectpuls((theta_norm - 1)/theta0);
Ry012 = (R0/(4*theta0^2)) .* rectpuls((theta_norm-0.5)/theta0);
Ry01 = Ry011 + Ry0111 + Ry012;

figure(2)
plot(theta_norm, Ry01)
title('Theoretical PSD')
xlabel('theta')
ylabel('Power Spectral Density')


%%

Rypm = (R0/(theta0)^2) .* rectpuls((theta_norm-0.5)/theta0);

figure(1)
plot(theta_norm, Rypm)
title('Theoretical PSD')
xlabel('theta')
ylabel('Power Spectral Density')
