%% Karan Mehta & Adnan Reddy, ECE5013 Project, 2022-04-11
% clc
% clear 
% close all

%% Physical Constants
c = 3e8;            % speed of light
T0 = 293;           % room temperature in Kelvin
K = 1.38e-23;       % Boltzman's constant

%% Parameters for MISO Radar System
fc = 10.5e9;            % center frequency = 10.5 GHz
lambda = c/fc;          % wavelength of radar system
BW = 125e6;             % total system bandwidth = 125 MHz
fp = 1e3;               % pulse repetition frequency = 1 kHz
Tp = 1/fp;              % pulse repetition interval = 1 ms
Np = 1;                % number of pulses


beta = 120e6;           % sweep bandwidth = 120 MHz
tau = 80e-6;            % pulse width = 80 usec

F = 10^(8/10);          % noise figure
Pt = 0.01*10^(20/10);   % transmit power
Gt = 10^(3/10);         % transmitter antenna gain
Gr = Gt;                % receiver antenna gain

theta_B = 60*pi/180;    % transmitter beamwidth

%% Parameters for Sampled System
fs = 250e6;         % sample rate = 250 Msamples/second
Ts = 1/fs;          % sample period

%% Time Vector
t = 0:Ts:Np*Tp-Ts;     % time vector (Np+1) Tp long, 

%% Parameters for Target
R0 = 30;            % R0 is initially 30 meters
theta = -10;        % theta = azimuth angle, initially at -10 degrees
v = 10;             % vertical  velocity is 10 m/s
RCS = 1;            % radar cross-section
%CPI = 881/64        % coherent processing interval

%% Received signals
za1 = zeros(size(t));   % received signal from transmitter 1
za2 = zeros(size(t));   % received signal from transmitter 2
za = zeros(size(t));

% Simulate noise
sigma_n = sqrt(K*T0*BW*F/2); % noise power
noise = sigma_n*(randn(1,length(t))) + 1i*(sigma_n*(randn(1,length(t))));

tx1=[ 0, lambda/4]; tx2=[0, -lambda/4]; rx=[0,0];

% calculate azimuth angle, range, and received signals at each location
for k=0:Np-1

    target=[R0*cosd(theta) R0*sind(theta)+k*Tp*v]; % default
%     target = [initial_x initial_y(CPI)+k*Tp*v]; % for CPI for loop, COMMENT OTHERWISE
    
    Rup1 = norm(tx1-target);    % Rup1 = distance between TX1 antenna and target
    Rup2 = norm(tx2-target);    % Rup2 = distance between TX2 antenna and target
    Rdown = norm(rx-target);    % Rdown = distance between RX antenna and target

    % za1 = received signal from TX1
    Ac = sqrt(Pt*Gt*Gr*lambda^2*RCS/(4*pi)^3/Rup1^2/Rdown^2);

    za1 = za1 + Ac*(rpulse(t-k*Tp-(Rup1+Rdown)/c,tau)) .* ...
        (exp(-1i*(2*pi/lambda)*(Rup1+Rdown))) .* ...
        (exp(1i*pi*(beta/tau).* (t-(tau/2)-(k*Tp)-(Rup1+Rdown)/c).^2  ));
    
    % za2 = received signal from TX2
    Ac = sqrt(Pt*Gt*Gr*lambda^2*RCS/(4*pi)^3/Rup2^2/Rdown^2);
    
    za2 = za2 + Ac*(rpulse(t-k*Tp-(Rup2+Rdown)/c,tau)) .* ...t
        (exp(-1i*(2*pi/lambda)*(Rup2+Rdown))) .* ...
        (exp(-1i*pi*(beta/tau).* (t-(tau/2)-(k*Tp)-(Rup2+Rdown)/c).^2  ));
    
    % za = total received signal
    za = za1 + za2;
end
za=za+noise;

%Prepare match filter for Tx1
tm=0:Ts:tau-Ts;
h1=exp(1i*pi*(beta/tau).*(tm-(tau/2)) .^2); 
h1=conj(fliplr(h1)); %conjugate and flip the time.

%Prepare match filter for Tx2
tm=0:Ts:tau-Ts;
h2=exp(-1i*pi*(beta/tau).*(tm-(tau/2)) .^2); 
h2=conj(fliplr(h2)); %conjugate and flip the time.

receivearray=zeros(Np,20500);
for k=0:Np-1
    receivearray(k+1,:)=za(k*uint64(Tp*fs)+1:round(k*uint64(Tp*fs)+tau*fs+500));
end

matcharray1=zeros(Np,501);
for k=0:Np-1
   matcharray1(k+1,:)=conv(receivearray(k+1,:),h1,'valid');
end

matcharray2=zeros(Np,501);
for k=0:Np-1
    matcharray2(k+1,:)=conv(receivearray(k+1,:),h2,'valid');
end

%%
M = 64;

taugrid=(0:500);
nugrid=1/M*(-M/2:1:(M/2)-1);

rangedoppler1=fftshift( fft(matcharray1,M,1),1);
rangedoppler2=fftshift( fft(matcharray2,M,1),1);

%% Graph

subplot(2,2,1)
imagesc(taugrid*Ts*c/2,1:Np,abs(matcharray1));
xlabel('Range (m)'); ylabel('Pulse No');

subplot(2,2,2)
imagesc(taugrid*Ts*c/2,nugrid*lambda*fp/2, abs(rangedoppler1))
xlabel('Range (m)'); ylabel('Velocity (m/sec)');

subplot(2,2,3)
imagesc(taugrid*Ts*c/2,1:Np,abs(matcharray2));
xlabel('Range (m)'); ylabel('Pulse No');

subplot(2,2,4)
imagesc(taugrid*Ts*c/2,nugrid*lambda*fp/2, abs(rangedoppler2))
xlabel('Range (m)'); ylabel('Velocity (m/sec)');

%Things to try, window the match filter, window the pulses, add ground
%return, use MTI cancelling, zeropad FFT, label the range doppler map in 
%velocity(m/sec) and Range( meters)

%% SNR
 
disp('theoretical')
disp(Pt*Gt*Gr*lambda^2*RCS/(4*pi)^3/Rdown^4/(K*T0*BW*F)) % estimate using Rdown to account for two different distances
disp('experimental')
disp(   (mean(abs(za1).^2)/mean(abs(noise).^2)) /   (fp*tau))
disp(   (mean(abs(za2).^2)/mean(abs(noise).^2)) /   (fp*tau))

%% Rectangular Pulse
function p = rpulse(t,tau)
    p = (t<=tau)&(t>=0);
end