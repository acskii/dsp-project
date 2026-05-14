% =========================================================================
% DSP FINAL PROJECT - PART II: Multi-Band Speech Equalizer
% -------------------------------------------------------------------------
% This script implements a graphic equalizer for speech enhancement.
% It allows the user to choose between a preset 7‑band mode and a custom
% mode with 5‑10 bands. For each band a digital filter (FIR Hamming or IIR
% Butterworth) is designed, the signal is filtered, a dB gain is applied,
% and all bands are summed. The result is compared to the original using
% time‑domain plots, magnitude spectra, Welch PSD, and spectrograms.
% The filtered audio is played and saved to a WAV file.
%
% The code uses only functions covered in the DSP lab materials plus basic
% MATLAB operations. Complex built‑in functions like pwelch or spectrogram
% are manually implemented to ensure full transparency.
% =========================================================================

clear; close all; clc;

%% ======================== USER INPUT SECTION ===========================
% In this section we ask the user for all necessary parameters.
% Every input is validated with while loops until correct values are given.

% ----- Audio file -----
fname = input('Enter the audio file name (e.g., ''speech.wav''): ', 's');
while ~exist(fname, 'file')
    fprintf('File "%s" not found. Please check the name and try again.\n', fname);
    fname = input('Enter the audio file name: ', 's');
end

% Read the audio file using audioread (modern replacement for wavread)
[x_orig, Fs_in] = audioread(fname);

% If stereo, convert to mono by averaging the two channels.
if size(x_orig,2) == 2
    x_orig = mean(x_orig, 2);
end

% ----- Output sample rate -----
Fs_out = input('Enter the desired output sample rate (Hz): ');
while ~isnumeric(Fs_out) || ~isscalar(Fs_out) || Fs_out <= 0
    fprintf('Invalid sample rate. It must be a positive number.\n');
    Fs_out = input('Enter the desired output sample rate (Hz): ');
end

% Resample the audio to the output sample rate if it differs from the input rate.
% Using rat to find a rational approximation of the ratio Fs_out/Fs_in.
if Fs_out ~= Fs_in
    tol = 1e-6;
    [P, Q] = rat(Fs_out / Fs_in, tol);
    x = resample(x_orig, P, Q);   % resample is covered in Lab 0 (DSP_Lab0)
    Fs = Fs_out;
else
    x = x_orig;
    Fs = Fs_in;
end
N_signal = length(x);
t = (0:N_signal-1)' / Fs;        % time vector for plotting

% ----- Operating mode -----
fprintf('\nSelect equalizer mode:\n');
fprintf('  1 - Preset (7 bands, speech-optimised)\n');
fprintf('  2 - Custom (5-10 bands, user-defined edges)\n');
mode = input('Enter 1 or 2: ');
while ~ismember(mode, [1,2])
    fprintf('Invalid mode. Please enter 1 or 2.\n');
    mode = input('Enter 1 or 2: ');
end

% Define band edges and gains
if mode == 1
    % Preset mode: 7 fixed bands (Hz)
    band_edges = [0, 100, 300, 800, 2000, 5000, 10000, 20000];
    num_bands = 7;
    fprintf('\nPreset bands (Hz):\n');
    for k = 1:num_bands
        fprintf('  Band %d: %d - %d\n', k, band_edges(k), band_edges(k+1));
    end
    % Ask for gains (one per band)
    gains_dB = input('Enter a vector of 7 gains (dB), e.g., [0,0,0,3,-2,0,0]: ');
    while ~isnumeric(gains_dB) || length(gains_dB) ~= 7
        fprintf('You must enter exactly 7 gain values.\n');
        gains_dB = input('Enter a vector of 7 gains (dB): ');
    end
else
    % Custom mode
    num_bands = input('Enter the number of bands (5-10): ');
    while ~isnumeric(num_bands) || ~isscalar(num_bands) || ...
          floor(num_bands) ~= num_bands || num_bands < 5 || num_bands > 10
        fprintf('Number of bands must be an integer between 5 and 10.\n');
        num_bands = input('Enter the number of bands (5-10): ');
    end
    % Band edges: the user must provide a vector of length num_bands+1,
    % starting with 0 and ending with 20000.
    fprintf('\nDefine band edges. The vector must start with 0 and end with 20000,\n');
    fprintf('and have exactly %d elements (e.g., [0 200 500 2000 8000 20000] for 5 bands).\n', num_bands+1);
    band_edges = input('Enter the band edges (Hz) as a row vector: ');
    while ~isnumeric(band_edges) || length(band_edges) ~= num_bands+1 || ...
          band_edges(1) ~= 0 || band_edges(end) ~= 20000 || ...
          any(diff(band_edges) <= 0)
        fprintf('Invalid band edges. Ensure:\n');
        fprintf('  - length is %d\n', num_bands+1);
        fprintf('  - first value is 0\n');
        fprintf('  - last value is 20000\n');
        fprintf('  - values are strictly increasing\n');
        band_edges = input('Enter the band edges (Hz) as a row vector: ');
    end
    % Gains
    gains_dB = input(sprintf('Enter a vector of %d gains (dB): ', num_bands));
    while ~isnumeric(gains_dB) || length(gains_dB) ~= num_bands
        fprintf('You must enter exactly %d gain values.\n', num_bands);
        gains_dB = input(sprintf('Enter a vector of %d gains (dB): ', num_bands));
    end
end

% ----- Filter type -----
fprintf('\nSelect filter type:\n');
fprintf('  FIR  - Finite Impulse Response (Hamming window)\n');
fprintf('  IIR  - Infinite Impulse Response (Butterworth)\n');
filter_type = upper(input('Enter FIR or IIR: ', 's'));
while ~ismember(filter_type, {'FIR','IIR'})
    fprintf('Invalid choice. Please enter FIR or IIR.\n');
    filter_type = upper(input('Enter FIR or IIR: ', 's'));
end

% ----- Filter order -----
order = input('Enter the filter order (positive integer): ');
while ~isnumeric(order) || ~isscalar(order) || order < 1 || floor(order) ~= order
    fprintf('Filter order must be a positive integer.\n');
    order = input('Enter the filter order: ');
end

%% =============== FILTER DESIGN FOR EACH BAND ===========================
% We treat each band as either a lowpass, bandpass, or highpass.
% The maximum frequency is limited by the Nyquist frequency (Fs/2).
% All band edges are clamped to be at most Fs/2.
Nyquist = Fs / 2;

% Prepare cell arrays to store filter coefficients
b_coeffs = cell(num_bands, 1);   % numerator coefficients
a_coeffs = cell(num_bands, 1);   % denominator coefficients (a = 1 for FIR)

fprintf('\n--- Filter Design Summary ---\n');
for idx = 1:num_bands
    f_low  = band_edges(idx);
    f_high = band_edges(idx+1);
    
    % Clamp frequencies to Nyquist
    f_low  = min(f_low, Nyquist);
    f_high = min(f_high, Nyquist);
    
    % Determine filter type and normalized cutoff(s)
    if f_low == 0
        % Lowpass
        ftype = 'low';
        Wn = f_high / Nyquist;
        descr = sprintf('Lowpass  0 - %.0f Hz', f_high);
    elseif f_high >= Nyquist
        % Highpass
        ftype = 'high';
        Wn = f_low / Nyquist;
        descr = sprintf('Highpass %.0f - %.0f Hz', f_low, Nyquist);
    else
        % Bandpass
        ftype = 'bandpass';
        Wn = [f_low, f_high] / Nyquist;
        descr = sprintf('Bandpass %.0f - %.0f Hz', f_low, f_high);
    end
    
    % Design the filter according to the chosen type
    if strcmp(filter_type, 'FIR')
        % FIR design using fir1 with Hamming window (Lab materials)
        win = hamming(order+1, 'periodic');   % hamming window of length order+1
        b = fir1(order, Wn, ftype, win);
        a = 1;
        fprintf('Band %d (%s): FIR order %d\n', idx, descr, order);
    else  % IIR Butterworth
        % Butterworth design. For bandpass/highpass, butter returns a filter
        % of order 2*order. This is standard behaviour (mentioned in labs).
        [b, a] = butter(order, Wn, ftype);
        fprintf('Band %d (%s): IIR Butterworth order %d\n', idx, descr, order);
    end
    
    b_coeffs{idx} = b;
    a_coeffs{idx} = a;
end

%% =============== FILTERING AND GAIN APPLICATION =========================
% We filter the (possibly resampled) input signal through each band filter,
% apply the gain, and sum all contributions.
y_bands = zeros(N_signal, num_bands);

for idx = 1:num_bands
    y_bands(:, idx) = filter(b_coeffs{idx}, a_coeffs{idx}, x);
    % Convert gain from dB to linear and apply
    gain_lin = 10^(gains_dB(idx)/20);
    y_bands(:, idx) = y_bands(:, idx) * gain_lin;
end

% Sum all bands to obtain the equalised signal
y_eq = sum(y_bands, 2);

% For FIR filters, the delay is exactly order/2 samples (constant group delay).
% We shift the output back to align with the original for fair comparison.
if strcmp(filter_type, 'FIR')
    delay_samples = floor(order/2);
    y_eq_aligned = y_eq(delay_samples+1 : end);
    x_aligned = x(1 : end-delay_samples);
    t_aligned = t(1 : end-delay_samples);
else
    % IIR: no simple alignment, we keep signals as is (will discuss phase distortion)
    y_eq_aligned = y_eq;
    x_aligned = x;
    t_aligned = t;
end

%% =============== ANALYSIS AND VISUALISATION =============================

% ----- 1. Time-domain comparison -----
figure('Name', 'Time-domain comparison');
plot(t_aligned, x_aligned, 'b'); hold on;
plot(t_aligned, y_eq_aligned, 'r');
xlabel('Time (s)'); ylabel('Amplitude');
title('Original (blue) vs. Equalised (red)');
legend('Original', 'Equalised');
grid on;

% ----- 2. Magnitude spectra (full signal) -----
NFFT = length(x_aligned);
X = abs(fftshift(fft(x_aligned)));   % magnitude spectrum of original
Y = abs(fftshift(fft(y_eq_aligned))); % magnitude spectrum of equalised
f_axis = linspace(-Fs/2, Fs/2, NFFT);

figure('Name', 'Magnitude Spectra');
plot(f_axis, 20*log10(X), 'b'); hold on;
plot(f_axis, 20*log10(Y), 'r');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Magnitude spectrum: original vs. equalised');
legend('Original', 'Equalised');
grid on; xlim([-Fs/2, Fs/2]);

% ----- 3. Welch Power Spectral Density (manual implementation) -----------
% We divide the signal into overlapping segments, apply a Hamming window,
% compute the periodogram of each segment, and average.
Nfft_welch = 256;          % FFT length
overlap   = 0.5;           % 50% overlap
window    = hamming(Nfft_welch, 'periodic');
% Compute PSD for original and equalised
[psd_orig, f_welch] = welch_manual(x_aligned, window, Nfft_welch, round(Nfft_welch*(1-overlap)), Fs);
[psd_eq,   ~]        = welch_manual(y_eq_aligned, window, Nfft_welch, round(Nfft_welch*(1-overlap)), Fs);

figure('Name', 'Welch PSD');
plot(f_welch, 10*log10(psd_orig), 'b'); hold on;
plot(f_welch, 10*log10(psd_eq), 'r');
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
title('Welch PSD: original (blue) vs. equalised (red)');
legend('Original', 'Equalised');
grid on;

% ----- 4. Spectrogram (manual STFT) --------------------------------------
% We compute short-time Fourier transform with overlapping windows.
win_len = 256;
overlap_stft = round(win_len * 0.75);
[S_orig, f_stft, t_stft] = stft_manual(x_aligned, win_len, overlap_stft, Fs);
[S_eq,   ~,       ~]      = stft_manual(y_eq_aligned, win_len, overlap_stft, Fs);

figure('Name', 'Spectrograms');
subplot(2,1,1);
imagesc(t_stft, f_stft, 20*log10(abs(S_orig))); axis xy; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Original spectrogram');
ylim([0, Fs/2]); caxis([-60 0]);   % set colour limits for clarity

subplot(2,1,2);
imagesc(t_stft, f_stft, 20*log10(abs(S_eq))); axis xy; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Equalised spectrogram');
ylim([0, Fs/2]); caxis([-60 0]);

% ----- 5. Filter characteristics for each band --------------------------
% For each band we plot: magnitude & phase response, impulse & step response,
% and pole-zero diagram.
for idx = 1:num_bands
    b = b_coeffs{idx};
    a = a_coeffs{idx};
    
    figure('Name', sprintf('Band %d characteristics', idx));
    
    % Magnitude and phase response (freqz)
    subplot(2,3,1);
    freqz(b, a, 1024, Fs);
    title(sprintf('Band %d: Magnitude & Phase', idx));
    
    % Impulse response
    subplot(2,3,2);
    impz(b, a, [], Fs);
    title('Impulse Response');
    
    % Step response (cumulative sum of impulse response)
    [h_imp, t_imp] = impz(b, a, [], Fs);
    step_resp = cumsum(h_imp);
    subplot(2,3,3);
    plot(t_imp, step_resp);
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Step Response'); grid on;
    
    % Pole-zero plot
    subplot(2,3,4);
    zplane(b, a);
    title('Pole-Zero Plot');
end

%% =============== PLAYBACK AND SAVE ======================================
% Play the equalised audio
fprintf('\nPlaying the equalised audio...\n');
sound(y_eq_aligned, Fs);   % ensure we use the aligned version (same length as original)
% Save the equalised signal to a WAV file
audiowrite('equalized_output.wav', y_eq_aligned, Fs);
fprintf('Equalised audio saved as "equalized_output.wav".\n');

%% =============== HELPER FUNCTIONS (manual implementations) ==============
% These are written from scratch to avoid using toolbox functions not explicitly
% covered in the labs.

function [psd, f] = welch_manual(x, window, nfft, noverlap, Fs)
% Manual Welch PSD estimate.
% x       - input signal (column vector)
% window  - analysis window (e.g., hamming(nfft,'periodic'))
% nfft    - FFT length
% noverlap- number of overlapping samples
% Fs      - sampling frequency
%
% The signal is divided into segments of length nfft, shifted by
% (nfft - noverlap) samples. Each segment is windowed, its periodogram
% computed, and the average is returned.
    win = window(:);   % ensure column
    seg_shift = nfft - noverlap;
    if seg_shift <= 0
        error('Overlap too large; segments must advance.');
    end
    % Number of segments
    n_seg = floor((length(x) - noverlap) / seg_shift);
    psd_sum = zeros(nfft, 1);
    for k = 1:n_seg
        start_idx = (k-1)*seg_shift + 1;
        seg = x(start_idx : start_idx + nfft - 1);
        seg_windowed = seg .* win;
        seg_fft = fft(seg_windowed, nfft);
        psd_sum = psd_sum + abs(seg_fft).^2;
    end
    psd = psd_sum / (n_seg * Fs * sum(win.^2));  % normalisation (Welch's method)
    psd = psd(1:floor(nfft/2)+1);                % keep positive frequencies
    psd(2:end-1) = 2*psd(2:end-1);               % compensate for single-sided
    f = linspace(0, Fs/2, length(psd));
end

function [S, f, t] = stft_manual(x, win_len, noverlap, Fs)
% Manual Short-Time Fourier Transform spectrogram.
% x        - input signal (column vector)
% win_len  - window length (samples)
% noverlap - number of overlapping samples
% Fs       - sampling frequency
%
% Returns:
% S  - complex STFT matrix (frequency x time)
% f  - frequency vector (Hz)
% t  - time vector (s)
    win = hamming(win_len, 'periodic');
    seg_shift = win_len - noverlap;
    if seg_shift <= 0
        error('Invalid overlap; seg_shift must be positive.');
    end
    n_seg = floor((length(x) - noverlap) / seg_shift);
    nfft = win_len;   % same length as window for simplicity
    S = zeros(nfft, n_seg);
    for k = 1:n_seg
        start_idx = (k-1)*seg_shift + 1;
        seg = x(start_idx : start_idx + win_len - 1);
        seg_windowed = seg .* win;
        S(:, k) = fft(seg_windowed, nfft);
    end
    f = linspace(0, Fs, nfft)';      % full frequency axis
    t = ((0:n_seg-1) * seg_shift + (win_len/2)) / Fs;  % time instants (centre of windows)
    % Keep only positive frequencies up to Nyquist for plotting
    keep = 1:floor(nfft/2)+1;
    S = S(keep, :);
    f = f(keep);
end