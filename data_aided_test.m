close all
clear variables

%% 1. System Parameters
on1 = 1;
fs = 10000;          % Sampling frequency (Hz)
T = 1/fs;           
numSymbols = 20;    % Increased for better lock visualization
sps = 8;            % Samples per symbol
symRate = fs/sps;    % Symbol rate (1000 Baud)
alpha = 0.5;         % Roll-off factor
filterLenSymbols = 10;
txFilter = comm.RaisedCosineTransmitFilter(...
    "FilterSpanInSymbols", filterLenSymbols, ...
    "RolloffFactor", alpha, "OutputSamplesPerSymbol", sps);
rxFilter = comm.RaisedCosineReceiveFilter(...
    "FilterSpanInSymbols", filterLenSymbols, ...
    "RolloffFactor", alpha, "InputSamplesPerSymbol", sps, ...
    "DecimationFactor", 1);

preamble = [-1 -1 -1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1]';
% Generate Random BPSK Data
data = [preamble; sign(randn(numSymbols, 1))];

if on1 == 0
    data = txFilter(data);
end

% Channel: Frequency and Phase Offset
t = (0:length(data)-1)'*T;
freq_offset = 20;    % 15 Hz offset
phase_offset = pi/4; % 45 degrees initial phase
input_signal = data .* exp(1i*(2*pi*freq_offset*t + phase_offset));

if on1 == 0
    input_signal = rxFilter(input_signal);
    SyncedData = dataAidedSync_not_symbol_sync(input_signal);
else

    SyncedData = dataAidedSync_symbol_sync(input_signal);
end


constDiagram = comm.ConstellationDiagram( ...
    SamplesPerSymbol=8, ...
    SymbolsToDisplaySource='Property', ...
    SymbolsToDisplay=100);
% Plotting the synchronized data  
%constDiagram(SyncedData)

figure;
subplot(2,1,1)
scatter(real(input_signal),imag(input_signal))
title("Input signal")
subplot(2,1,2)
scatter(real(SyncedData), imag(SyncedData))
title("Synchronized signal")
%dataAidedSync(testSym.')

function syncdSymOut = dataAidedSync_symbol_sync(symIn)
%symRate = 1000;
%sample_rate = 10;
%global T;

% symIn: preamble sequence + data symbols
% For use in Simulink model 
persistent cfo_est 
alpha = 0.7;

if isempty(cfo_est)
    cfo_est = 0;
end

preamble = [-1 -1 -1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1].';

if coder.target('MATLAB')
    % Standard MATLAB logic
    r_preamb = symIn((1:length(preamble)));
else
    % Logic specific to Simulink Execution
    r_preamb = symIn((1:length(preamble)+10));
end

% Use the functiomn movemean() to lowpass filter the extracted phase

L = length(preamble);
z = r_preamb.* conj(preamble);
M = ceil(L/1.5);
R = zeros(M,1);
for k=1:M
    R(k) = sum(z(k+1:end) .* conj(z(1:end-k)));
end

dw = mean(angle(R) ./ (1:M)')

t = (0:L-1)';
z_corrected = z .* exp(-1j * dw * t);
phi0 = angle(sum(z_corrected))


% For Simulink model use: 
cfo_est = dw*alpha + (1-alpha)*cfo_est;

%cfo_est = 10;

% Rotate Data to compensate for phase and frequency offset
% NB! For use in the simulink model (again) remember that there is a sample
% offset of 10 samples before the preamble starts (due to application
% specific reasons) 
n = (L:length(symIn)-1).';
syncdSymOut = symIn(L:length(symIn)-1) .* exp(1j*(-dw*n - phi0));

end

function syncdSymOut = dataAidedSync_not_symbol_sync(symIn)
%symRate = 1000;
%sample_rate = 10;
%global T;

% symIn: preamble sequence + data symbols
% For use in Simulink model 
persistent cfo_est 
alpha = 0.7;

if isempty(cfo_est)
    cfo_est = 0;
end
filterLenSymbols = 10;
sps = 8;

txFilter = comm.RaisedCosineTransmitFilter(...
    "FilterSpanInSymbols", filterLenSymbols, ...
    "RolloffFactor", alpha, "OutputSamplesPerSymbol", sps, "Shape","Normal");
preamble = txFilter([-1 -1 -1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1].');

if coder.target('MATLAB')
    % Standard MATLAB logic
    r_preamb = symIn((1:length(preamble)));
else
    % Logic specific to Simulink Execution
    r_preamb = symIn((1:length(preamble)+10));
end

% Use the functiomn movemean() to lowpass filter the extracted phase

L = length(preamble);
z = r_preamb.* conj(preamble);
M = ceil(L/2);
R = zeros(M,1);
for k=1:M
    R(k) = sum(z(k+1:end) .* conj(z(1:end-k)));
end

dw = mean(angle(R) ./ (1:M)')

t = (0:L-1)';
z_corrected = z .* exp(-1j * dw * t);
phi0 = angle(sum(z_corrected))


% For Simulink model use: 
cfo_est = dw*alpha + (1-alpha)*cfo_est;

%cfo_est = 10;

% Rotate Data to compensate for phase and frequency offset
% NB! For use in the simulink model (again) remember that there is a sample
% offset of 10 samples before the preamble starts (due to application
% specific reasons) 
n = (L:length(symIn)-1).';
syncdSymOut = symIn(L:length(symIn)-1) .* exp(1j*(-dw*n - phi0));

end




