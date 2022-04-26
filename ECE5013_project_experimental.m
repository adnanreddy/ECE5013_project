%% Karan Mehta & Adnan Reddy, ECE5013 Project, 2022-04-25
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
Np = 64;                % number of pulses


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

%%

M=256; % number of pulses in one CPI
N_cpi= 32; % number of CPIs in data

figure(800);
while size(findobj(800)) > 0
    figure(800);
    fid=fopen('measured_data.bin','rb');
    for ii=1:N_cpi
    % for ii=N_cpi:N_cpi % for testing
        phase_history1=zeros(301,M);
        phase_history2=zeros(301,M);

        for jj=1:M
            phase_history1(:,jj)=fread(fid,301,'double');
            phase_history1(:,jj)=phase_history1(:,jj)+1i*fread(fid,301,'double');

            phase_history2(:,jj)=fread(fid,301,'double');
            phase_history2(:,jj)=phase_history2(:,jj)+1i*fread(fid,301,'double');
        end

        %only Range bins >=131 have target information, earlier range bins are
        %due to antenna coupling, TX leaking into RX, delays in the radar.
        % row range-bins 0:170, columns :pulse_no 1:M

        phase_history1=phase_history1(131:end,:);
        phase_history2=phase_history2(131:end,:);

        phase_history1 = phase_history1';
        phase_history2 = phase_history2';

        % phase_history matrix: row range-bins 0:170, columns :pulse_no 1:M
        % Construct Range Doppler map from the two antennas and display.
        % Consider canceling stationary clutter (two/three pulse canceller, or 
        % subtracting average (across slow time) range-profile from each return

        % SUBTRACT SLOW TIME AVERAGE
        phase_history1=phase_history1 - mean(phase_history1,1);
        phase_history2=phase_history2 - mean(phase_history2,1);

        % TWO PULSE CANCELLER

%         phase_history1 = phase_history1 - [zeros()]

        % GRAPH
        taugrid=(130:300);
        nugrid=1/M*(-M/2:1:(M/2)-1);


        rangedoppler1=fftshift( fft(hamming(M).*phase_history1,M,1),1);
        rangedoppler2=fftshift( fft(hamming(M).*phase_history2,M,1),1);

        subplot(2,2,1)
        imagesc(taugrid*Ts*c/2,1:Np,abs(phase_history1));
        xlabel('Range (m)'); ylabel('Pulse No');
        title_string = strcat('Slow-Fast Time, Tx1, CPI=',num2str(ii)); title(title_string);

        subplot(2,2,3)
        imagesc(taugrid*Ts*c/2,nugrid*lambda*fp/2, abs(rangedoppler1))
        xlabel('Range (m)'); ylabel('Velocity (m/sec)');
        title_string = strcat('Range-Doppler, Tx1, CPI=',num2str(ii)); title(title_string);

        subplot(2,2,2)
        imagesc(taugrid*Ts*c/2,1:Np,abs(phase_history2));
        xlabel('Range (m)'); ylabel('Pulse No');
        title_string = strcat('Slow-Fast Time, Tx2, CPI=',num2str(ii)); title(title_string);

        subplot(2,2,4)
        imagesc(taugrid*Ts*c/2,nugrid*lambda*fp/2, abs(rangedoppler2))
        xlabel('Range (m)'); ylabel('Velocity (m/sec)');
        title_string = strcat('Range-Doppler, Tx2, CPI=',num2str(ii)); title(title_string);
        
        %% Range and Velocity

        ard1 = abs(rangedoppler1); ard2 = abs(rangedoppler2); % absolute range doppler
        [rows,cols] = size(ard1); % if they aren't the same size I'm out of luck

        peaks = [0 0];
        index = [1 1 1 1];

        for row = 1:rows
            for col = 1:cols
                if peaks(1) < ard1(row,col)
                    peaks(1) = ard1(row,col);
                    index(1:2) = [row,col];
                end
                if peaks(2) < ard2(row,col)
                    peaks(2) = ard2(row,col);
                    index(3:4) = [row,col];
                end
            end
        end

        range_CPI_1(ii) = taugrid(index(2))*Ts*c/2;
        velocity_CPI_1(ii) = nugrid(index(1))*lambda*fp/2;
        range_CPI_2(ii) = taugrid(index(4))*Ts*c/2;
        velocity_CPI_2(ii) = nugrid(index(3))*lambda*fp/2;
        
        % Angle

        for jj = 1:Np
            [peak,index1] = max(abs(matcharray1(jj,:))); % find the range bin with highest magnitude return
            [peak,index2] = max(abs(matcharray2(jj,:)));

            angle_pulse(jj) = (angle(matcharray1(jj,index1) / matcharray2(jj,index2))); % absolute phase difference
        end
        angle_CPI(ii) = mean(asind(angle_pulse)*(lambda/2/pi)/(lambda/2)); % see report for diagram
        
        pause(0.5);

    end
    fclose(fid);
    
    figure(801); plot(range_CPI_1); hold on; plot(range_CPI_2); title('Range'); xlabel('CPI'); ylabel('Meters'); legend('Tx1','Tx2'); hold off; saveas(gcf,'range_exp.jpg')
    figure(802); plot(velocity_CPI_1); hold on; plot(velocity_CPI_2); title('Velocity'); xlabel('CPI'); ylabel('Meters per Second'); hold off; saveas(gcf,'velocity_exp.jpg')
    figure(803); scatter(1:N_cpi,angle_CPI); title('Azimuth Angle'); xlabel('CPI'); ylabel('Degrees'); saveas(gcf,'angle_exp.jpg')
end