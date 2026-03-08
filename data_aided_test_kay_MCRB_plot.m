clear all
close all

fs = 1000;          % Sampling frequency (Hz)
T = 1/fs;           
numSymbols = 100;    % Increased for better lock visualization
sps = 1;            % Samples per symbol
symRate = fs/sps;    % Symbol rate (1000 Baud)

preamble = [-1 -1 -1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1]';

t = (0:numSymbols + length(preamble)-1)'*T;
freq_offset = 15;    % 15 Hz offset
phase_offset = pi/4; % 45 degrees initial phase
N = 300;
phi0 = zeros(N, 1);
dw = zeros(N,1);
SNR_vec = 0:0.5:30;
SNR_lin_vec = 10.^(SNR_vec./10);
variance_f = zeros(length(SNR_vec),1);
variance_phi = zeros(length(SNR_vec),1);
MCRB_f = (3/(8*(pi^2)*(length(preamble)^3)*T^2))./(SNR_lin_vec);
MCRB_phi = (1/(2*(length(preamble))))./(SNR_lin_vec);
for j = 1:length(SNR_lin_vec)
for i = 1:N
    data = [preamble; sign(randn(numSymbols, 1))];
    input_signal = data .* exp(1i*(2*pi*freq_offset*t + phase_offset));
    w_offset = 2*pi*freq_offset;
    input_signal = awgn(input_signal, SNR_vec(j));
    [SyncedData, dw(i), phi0(i)] = dataAidedSync_symbol_sync(input_signal);
    
end
    variance_f(j) = var(dw*fs/(2*pi));
    variance_phi(j) = var(phi0);

end

figure
subplot(2,1,1)
semilogy(SNR_vec, variance_f, LineWidth=1.5, DisplayName="Computed variance")
hold on
grid on
semilogy(SNR_vec, MCRB_f, LineWidth=1.5, DisplayName="MCRB")
legend
title("Frequency variance")
xlabel("SNR [dB]")
ylabel("Var[\Delta f]")
ax = gca();
ax.FontSize = 20;
subplot(2,1,2)

semilogy(SNR_vec, variance_phi, LineWidth=1.5, DisplayName="Computed variance")
hold on
semilogy(SNR_vec, MCRB_phi, LineWidth=1.5, DisplayName="MCRB")
grid on
legend
xlabel("SNR [dB]")
ylabel("Var[\Delta \phi]")
title("Phase variance")
ax = gca();
ax.FontSize = 20;


function [syncdSymOut, corr_freq , phi0] = dataAidedSync_symbol_sync(symIn) 
% Data aided sync using samples after symbol synchronization
% Optimized for a single 13-bit Barker sequence using Kay's Estimator

persistent cfo_est 
% Lower alpha means more smoothing. If CFO changes rapidly, increase this.
% If CFO is relatively stable, a smaller alpha (e.g., 0.1 to 0.3) helps hide the 13-bit noise.
alpha = 0.5; 

% Single 13-bit Barker sequence
preamble = [-1 -1 -1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1].';
L = length(preamble);

% --- 1. PREAMBLE ALIGNMENT & EXTRACTION ---
% Ensure we account for the 10-sample offset explicitly to align arrays perfectly.
if coder.target('MATLAB')
    offset = 0;
else
    offset = 10; 
end

% Protect against array bounds issues in Simulink
if length(symIn) <= offset + L
    syncdSymOut = symIn;
    return;
end

% Extract EXACTLY L samples perfectly aligned with the preamble
r_preamb = symIn(offset + 1 : offset + L);

% Remove the BPSK modulation
z = r_preamb .* conj(preamble);

% --- 2. IMPROVED CFO ESTIMATION (Kay's Estimator) ---
% Delay and multiply by 1 sample
diff_z = z(2:end) .* conj(z(1:end-1));
angles = angle(diff_z);

% Compute optimal parabolic weights for L samples
k = (1:L-1).';
weights = (1.5 * L / (L^2 - 1)) .* (1 - ((2*k - L) ./ L).^2);

% Instantaneous CFO estimate (radians/sample)
dw = sum(weights .* angles);

% --- 3. APPLY SMOOTHING ---
% Update the persistent variable
if isempty(cfo_est)
    cfo_est = dw;
else
    cfo_est = dw * alpha + (1 - alpha) * cfo_est;
end

% --- 4. PHASE ESTIMATION ---
% Remove the estimated frequency offset from the preamble to find the static phase offset.
% IMPORTANT: Use the smoothed cfo_est here, not the instantaneous dw.
t_preamb = (0:L-1).';
z_corrected = z .* exp(-1j * cfo_est * t_preamb);
phi0 = mean(angle(z_corrected));

corr_freq = cfo_est;
% --- 5. DATA CORRECTION ---
% Define indices for the data portion strictly AFTER the preamble
idx_data = (offset + L + 1 : length(symIn)).';

% Time vector for the data symbols relative to the start of the preamble
t_data = (L : length(symIn) - offset - 1).'; 

% Extract and rotate the data
data_in = symIn(idx_data);
syncdSymOut = data_in .* exp(-1j * (cfo_est * t_data + phi0));

end
