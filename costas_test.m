clear all; close all;

%% 1. System Parameters
fs = 10000;          % Sampling frequency (Hz)
T = 1/fs;           
numSymbols = 500;    % Increased for better lock visualization
sps = 10;            % Samples per symbol
symRate = fs/sps;    % Symbol rate (1000 Baud)
alpha = 0.5;         % Roll-off factor
filterLenSymbols = 8;

% Generate Random BPSK Data
data = sign(randn(numSymbols, 1));

% RRC Transmit Filter
txFilter = comm.RaisedCosineTransmitFilter(...
    "FilterSpanInSymbols", filterLenSymbols, ...
    "RolloffFactor", alpha, "OutputSamplesPerSymbol", sps);
data_upsampled = txFilter(data);

% Channel: Frequency and Phase Offset
t = (0:length(data_upsampled)-1)'*T;
freq_offset = 15;    % 15 Hz offset
phase_offset = pi/4; % 45 degrees initial phase
input_signal = data_upsampled .* exp(1i*(2*pi*freq_offset*t + phase_offset));

%% 2. Receiver: Matched Filtering (RRC Receive Filter)
% It is crucial to filter BEFORE the loop to minimize ISI
rxFilter = comm.RaisedCosineReceiveFilter(...
    "FilterSpanInSymbols", filterLenSymbols, ...
    "RolloffFactor", alpha, "InputSamplesPerSymbol", sps, ...
    "DecimationFactor", 1);
received_filtered = rxFilter(input_signal);

%% 3. Costas Loop Implementation
% Loop Design Parameters
loop_bw = 0.02 * symRate; % Normalized Loop Bandwidth
zeta = 0.707;             % Damping factor (critically damped)
kpd = 1;                  % Phase Detector Gain (assuming AGC is used)
k0 = 1;                   % VCO Gain

% Calculating PI Controller Gains
wn = (4 * loop_bw * zeta) / (zeta + 1/(4*zeta));
kp = (2 * zeta * wn * T) / (kpd * k0);
ki = (wn^2 * T^2) / (kpd * k0);

% Initialize Buffers
vco_phase = 0;
loop_integral = 0;
recovered_I = zeros(size(received_filtered));
error_track = zeros(size(received_filtered));

% Processing Loop
for n = 1:length(received_filtered)
    % Phase Correction (Complex Mixing)
    % We rotate the signal back by the estimated VCO phase
    sample = received_filtered(n) * exp(-1i * vco_phase);
    
    I = real(sample);
    Q = imag(sample);
    
    % Decision-Directed Phase Error Detector (Optimized for RRC)
    % Using sign(I) makes it robust against the amplitude ripples of RRC pulses
    error = sign(I) * Q; 
    
    % Loop Filter (PI Controller)
    loop_integral = loop_integral + ki * error;
    vco_speed = kp * error + loop_integral;
    
    % Update VCO Phase for the next sample
    vco_phase = vco_phase + vco_speed;
    
    % Store outputs
    recovered_I(n) = I;
    error_track(n) = error;
end

%% 4. Visualization
figure('Name', 'Costas Loop Performance');
subplot(2,1,1); 
plot(error_track); 
title('Phase Error (Discriminator Output)'); 
xlabel('Samples'); ylabel('Error'); grid on;

subplot(2,1,2); 
plot(recovered_I); 
title('Recovered In-Phase Signal (Baseband)'); 
xlabel('Samples'); ylabel('Amplitude'); grid on;

% Constellation Diagram
cd = comm.ConstellationDiagram('Title', 'Locked Constellation (After Transients)');
% Display only after the loop has likely reached a lock (skip first 2000 samples)
cd(recovered_I(2000:end)); 

figure('Name','Input compare');
plot(recovered_I)
hold on
plot(real(rxFilter(data_upsampled))) 
title('Recovered In-Phase Signal (Baseband)'); 
xlabel('Samples'); ylabel('Amplitude'); grid on;
