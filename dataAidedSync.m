function syncdSymOut = dataAidedSync(symIn)

% Data aided sync using samples after symbol synchronization
% Optimized for a single 13-bit Barker sequence using Kay's Estimator

persistent cfo_est 
% Lower alpha means more smoothing. If CFO changes rapidly, increase this.
% If CFO is relatively stable, a smaller alpha (e.g., 0.1 to 0.3) helps hide the 13-bit noise.
alpha = 1; 

if isempty(cfo_est)
    cfo_est = 0;
end

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
cfo_est = dw * alpha + (1 - alpha) * cfo_est;

% --- 4. PHASE ESTIMATION ---
% Remove the estimated frequency offset from the preamble to find the static phase offset.
% IMPORTANT: Use the smoothed cfo_est here, not the instantaneous dw.
t_preamb = (0:L-1).';
z_corrected = z .* exp(-1j * cfo_est * t_preamb);
phi0 = angle(sum(z_corrected));

% --- 5. DATA CORRECTION ---
% Define indices for the data portion strictly AFTER the preamble
idx_data = (offset + L + 1 : length(symIn)).';

% Time vector for the data symbols relative to the start of the preamble
t_data = (L : length(symIn) - offset - 1).'; 

% Extract and rotate the data
data_in = symIn(idx_data);
syncdSymOut = data_in .* exp(-1j * (cfo_est * t_data + phi0));


end