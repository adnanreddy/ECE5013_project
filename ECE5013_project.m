%% Karan Mehta & Adnan Reddy, ECE5013 Project, 2022-04-11
clc
clear 
close all

%% Physical Constants
c = 3e8;            % speed of light
T0 = 293;           % room temperature in Kelvin
k = 1.38e-23;       % Boltzman's constant
F = 10^(8/10);      % noise figure

%% Parameters for MISO Radar System
fc = 10.5e9;        % center frequency = 10.5 GHz
lambda = c/fc;      % wavelength of radar system
BW = 125e6;         % total system bandwidth = 125 MHz
fp = 1e3;           % pulse repetition frequency = 1 kHz
Tp = 1/fp;          % pulse repetition interval = 1 ms
Np = 64;            % number of pulses


beta = 120e6;       % sweep bandwidth = 120 MHz
tau = 80e-6;        % pulse width = 80 usec



% Here I assume SNR=0 dB, you have to simulate correct Amplitude and Noise
% Levels, using radar range equation, etc


%% Parameters for Sampled System
fs = 250e6;         % sample rate = 250 Msamples/second
Ts = 1/fs;          % sample period

%% Time Vector

t = 0:Ts:Np*Tp-Ts;     % time vector (Np+1) Tp long, 

%% Parameters for Target
R0 = 30;            % R0 is initially 30 meters
theta = -10;        % theta = azimuth angle, initially at -10 degrees
v = 10;             % vertical  velocity is 10 m/s

%% Received signals
za1 = zeros(size(t));   % received signal from transmitter 1
za2 = zeros(size(t));   % received signal from transmitter 2
za = zeros(size(t));

% Simulate noise
sigma_n = k*T0*BW*F; % noise power
noise = sigma_n/2*(randn(1,length(t))) + 1i*(sigma_n/2*(randn(1,length(t))));

tx1=[ 0, lambda/4]; tx2=[0, -lambda/4]; rx=[0,0];

% calculate azimuth angle, range, and received signals at each location
for k=0:Np-1

    target=[R0*cosd(theta) R0*sind(theta)]; %TODO add back velocity
    Rup1 = norm(tx1-target);    % Rup1 = distance between TX1 antenna and target
    Rup2 = norm(tx2-target);    % Rup2 = distance between TX2 antenna and target
    Rdown = norm(rx-target);    % Rdown = distance between RX antenna and target

	Ac = 1; % TODO need RCS for range equation?

    za1 = za1 + Ac*(rpulse(t-k*Tp-(Rup1+Rdown)/c,tau)) .* ...
        (exp(-1i*(2*pi/lambda)*(Rup1+Rdown))) .* ...
        (exp(1i*pi*(beta/tau).* (t-(tau/2)-(k*Tp)-(Rup1+Rdown)/c).^2  ));
    
    za2 = za2 + Ac*(rpulse(t-k*Tp-(Rup2+Rdown)/c,tau)) .* ...
        (exp(-1i*(2*pi/lambda)*(Rup2+Rdown))) .* ...
        (exp(1i*pi*(beta/tau).* (t-(tau/2)-(k*Tp)-(Rup2+Rdown)/c).^2  ));

    za = za1 + za2;
    
    % za1 = received signal from TX1
    % za2 = received signal from TX2
    % za = total received signal
end
za=za+noise;

%Prepare match filter for Tx1
tm=0:Ts:tau-Ts;
h1=exp(1i*pi*(beta/tau).*(tm-(tau/2)) .^2); 
h1=conj(fliplr(h1)); %conjugate and flip the time.


receivearray=zeros(Np,20500);
for k=0:Np-1
    receivearray(k+1,:)=za(k*uint64(Tp*fs)+1:round(k*uint64(Tp*fs)+tau*fs+500));
end


matcharray=zeros(Np,501);
for k=0:Np-1
    matcharray(k+1,:)=conv(receivearray(k+1,:),h1,'valid');
end

taugrid=(0:500)*Ts;
figure(1);imagesc(taugrid,1:64,abs(matcharray));
xlabel('Delay'); ylabel('Pulse No');

rangedoppler=fftshift( fft(matcharray,[],1),1);

nugrid=1/Np*(-Np/2:1:(Np/2)-1);
taugrid=(0:500)*Ts;
figure(2);imagesc(taugrid,nugrid, abs(rangedoppler))
xlabel('Delay (sec)'); ylabel('Normalized Frequency (sec)');

%Things to try, window the match filter, window the pulses, add ground
%return, use MTI cancelling, zeropad FFT, label the range doppler map in 
%velocity(m/sec) and Range( meters)


%% Rectangular Pulse
function p = rpulse(t,tau)
    p = (t<=tau)&(t>=0);
end