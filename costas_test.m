clear all; close all;

%% 1. System Parameters
fs = 10000;          % Sampling frequency (Hz)
T = 1/fs;           
numSymbols = 100;    % Increased for better lock visualization
sps = 8;            % Samples per symbol
symRate = fs/sps;    % Symbol rate (1000 Baud)
alpha = 0.5;         % Roll-off factor
fls = 10;             % Filter length in symbols

% Generate Random BPSK Data
data = [sign(randn(numSymbols -13, 1)) ];

% RRC Transmit Filter
txFilter = comm.RaisedCosineTransmitFilter(...
    "FilterSpanInSymbols", fls, ...
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
    "FilterSpanInSymbols", fls, ...
    "RolloffFactor", alpha, "InputSamplesPerSymbol", sps, ...
    "DecimationFactor", 1);
received_filtered = rxFilter(input_signal);
received_filtered = received_filtered(2*fls+1:end);
%% 3. Costas Loop Implementation
% Loop Design Parameters
loop_bw = 0.02 * symRate; % Normalized Loop Bandwidth
zeta = 0.707;             % Damping factor (critically damped)
kpd = 1;                  % Phase Detector Gain (assuming AGC is used)
k0 = 1;                   % VCO/NCO Gain

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
plot(error_track, "LineWidth",1.5); 
title('Phase Error (Discriminator Output)'); 
xlabel('Samples'); ylabel('Error'); grid on;
ax = gca();
ax.FontSize = 20;

% Constellation Diagram
cd = comm.ConstellationDiagram('Title', 'Locked Constellation (After Transients)');
% Display only after the loop has likely reached a lock (skip first 2000 samples)
cd(recovered_I(2000:end)); 

figure('Name','Input compare');
plot(recovered_I,"DisplayName", "Compensated samples", "LineWidth",1.5)
hold on
wo_freq_off = real(rxFilter(data_upsampled));
plot(wo_freq_off(2*fls+1:end),"DisplayName", "Transmitted samples wo frequency offset", "LineWidth",1.5) 
title('Recovered In-Phase Signal (Baseband)'); 
xlabel('Samples'); ylabel('Amplitude'); grid on;
ax = gca();
ax.FontSize = 20;
legend 
