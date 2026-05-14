% =========================================================================
% DSP FINAL PROJECT - PART II: Multi-Band Speech Equalizer
% -------------------------------------------------------------------------
% This script implements a graphic equalizer for speech enhancement.
% It allows the user to choose between a preset 7-band mode and a custom
% mode with 5-10 bands. For each band a digital filter (FIR Hamming or IIR
% Butterworth) is designed, the signal is filtered, a dB gain is applied,
% and all bands are summed. The result is compared to the original using
% time-domain plots, magnitude spectra, Welch PSD, and spectrograms.
% The filtered audio is played and saved to a WAV file.
%
% The code uses only functions covered in the DSP lab materials plus basic
% MATLAB operations. Complex built-in functions like pwelch or spectrogram
% are manually implemented to ensure full transparency.
% =========================================================================

% FIX 1: Removed duplicate 'clear' statement from original code.
clear; close all; clc;

%% ======================== USER INPUT SECTION ===========================

% ----- Audio file -----
fname = input('Enter the audio file name (e.g., ''speech.wav''): ', 's');
while ~exist(fname, 'file')
    fprintf('File "%s" not found. Please check the name and try again.\n', fname);
    fname = input('Enter the audio file name: ', 's');
end

[x_orig, Fs_in] = audioread(fname);

% If stereo, convert to mono by averaging the two channels.
if size(x_orig, 2) == 2
    x_orig = mean(x_orig, 2);
end

% ----- Output sample rate -----
Fs_out = input('Enter the desired output sample rate (Hz): ');
while ~isnumeric(Fs_out) || ~isscalar(Fs_out) || Fs_out <= 0
    fprintf('Invalid sample rate. It must be a positive number.\n');
    Fs_out = input('Enter the desired output sample rate (Hz): ');
end

% Resample if needed using rational approximation.
if Fs_out ~= Fs_in
    tol = 1e-6;
    [P, Q] = rat(Fs_out / Fs_in, tol);
    x = resample(x_orig, P, Q);
    Fs = Fs_out;
else
    x = x_orig;
    Fs = Fs_in;
end
N_signal = length(x);
t = (0:N_signal-1)' / Fs;

% ----- Operating mode -----
fprintf('\nSelect equalizer mode:\n');
fprintf('  1 - Preset (7 bands, speech-optimised)\n');
fprintf('  2 - Custom (5-10 bands, user-defined edges)\n');
mode = input('Enter 1 or 2: ');
while ~ismember(mode, [1, 2])
    fprintf('Invalid mode. Please enter 1 or 2.\n');
    mode = input('Enter 1 or 2: ');
end

if mode == 1
    band_edges = [0, 100, 300, 800, 2000, 5000, 10000, 20000];
    num_bands = 7;
    fprintf('\nPreset bands (Hz):\n');
    for k = 1:num_bands
        fprintf('  Band %d: %d - %d\n', k, band_edges(k), band_edges(k+1));
    end
    gains_dB = input('Enter a vector of 7 gains (dB), e.g., [0,0,0,3,-2,0,0]: ');
    while ~isnumeric(gains_dB) || length(gains_dB) ~= 7
        fprintf('You must enter exactly 7 gain values.\n');
        gains_dB = input('Enter a vector of 7 gains (dB): ');
    end
else
    num_bands = input('Enter the number of bands (5-10): ');
    while ~isnumeric(num_bands) || ~isscalar(num_bands) || ...
          floor(num_bands) ~= num_bands || num_bands < 5 || num_bands > 10
        fprintf('Number of bands must be an integer between 5 and 10.\n');
        num_bands = input('Enter the number of bands (5-10): ');
    end
    fprintf('\nDefine band edges. The vector must start with 0 and end with 20000,\n');
    fprintf('and have exactly %d elements.\n', num_bands+1);
    band_edges = input('Enter the band edges (Hz) as a row vector: ');
    while ~isnumeric(band_edges) || length(band_edges) ~= num_bands+1 || ...
          band_edges(1) ~= 0 || band_edges(end) ~= 20000 || ...
          any(diff(band_edges) <= 0)
        fprintf('Invalid band edges. Ensure length=%d, starts at 0, ends at 20000, strictly increasing.\n', num_bands+1);
        band_edges = input('Enter the band edges (Hz) as a row vector: ');
    end
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
while ~ismember(filter_type, {'FIR', 'IIR'})
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
Nyquist = Fs / 2;

% Wn must stay strictly inside (0, 1) for fir1 and butter.
% We use a small epsilon to keep values off the boundaries.
WN_EPS = 1e-6;

b_coeffs   = cell(num_bands, 1);
a_coeffs   = cell(num_bands, 1);
band_active = true(num_bands, 1);   % tracks which bands are actually designed

fprintf('\n--- Filter Design Summary ---\n');
for idx = 1:num_bands
    f_low  = band_edges(idx);
    f_high = band_edges(idx+1);

    % Clamp both edges to the Nyquist frequency.
    f_low  = min(f_low,  Nyquist);
    f_high = min(f_high, Nyquist);

    % FIX 9: Skip bands that become degenerate after clamping.
    % This happens when the sample rate is lower than the band edges
    % (e.g., Fs=400 Hz, Nyquist=200 Hz collapses any band above 200 Hz).
    if f_low >= f_high
        fprintf('Band %d (%.0f - %.0f Hz): SKIPPED (above Nyquist = %.0f Hz)\n', ...
                idx, band_edges(idx), band_edges(idx+1), Nyquist);
        % Store a pass-through (gain-only) filter so the indexing stays consistent.
        b_coeffs{idx}  = 1;
        a_coeffs{idx}  = 1;
        band_active(idx) = false;
        continue;
    end

    if f_low == 0
        ftype = 'low';
        Wn    = f_high / Nyquist;
        descr = sprintf('Lowpass  0 - %.0f Hz', f_high);
    elseif f_high >= Nyquist
        ftype = 'high';
        Wn    = f_low / Nyquist;
        descr = sprintf('Highpass %.0f - %.0f Hz', f_low, Nyquist);
    else
        ftype = 'bandpass';
        Wn    = [f_low, f_high] / Nyquist;
        descr = sprintf('Bandpass %.0f - %.0f Hz', f_low, f_high);
    end

    % FIX 10: Clamp Wn strictly away from 0 and 1.
    % fir1 and butter both reject Wn = 0 or Wn = 1 exactly.
    if isscalar(Wn)
        Wn = max(WN_EPS, min(1 - WN_EPS, Wn));
    else
        Wn(1) = max(WN_EPS,     Wn(1));
        Wn(2) = min(1 - WN_EPS, Wn(2));
        % If clamping made the two edges identical, skip this band too.
        if Wn(1) >= Wn(2)
            fprintf('Band %d (%s): SKIPPED (band too narrow after Nyquist clamping)\n', idx, descr);
            b_coeffs{idx}  = 1;
            a_coeffs{idx}  = 1;
            band_active(idx) = false;
            continue;
        end
    end

    if strcmp(filter_type, 'FIR')
        fir_order = order;
        % FIX 3: fir1 requires an EVEN order for bandpass filters.
        if strcmp(ftype, 'bandpass') && mod(fir_order, 2) ~= 0
            fir_order = fir_order + 1;
            fprintf('  [Note] Band %d: FIR order incremented to %d (must be even for bandpass).\n', idx, fir_order);
        end
        win = hamming(fir_order + 1);
        b   = fir1(fir_order, Wn, ftype, win);
        a   = 1;
        fprintf('Band %d (%s): FIR order %d\n', idx, descr, fir_order);
    else
        [b, a] = butter(order, Wn, ftype);
        fprintf('Band %d (%s): IIR Butterworth order %d\n', idx, descr, order);
    end

    b_coeffs{idx} = b;
    a_coeffs{idx} = a;
end

%% =============== FILTERING AND GAIN APPLICATION =========================
y_bands = zeros(N_signal, num_bands);

for idx = 1:num_bands
    % FIX 9 continued: do not filter skipped (degenerate) bands —
    % their contribution stays zero so they don't corrupt the output.
    if ~band_active(idx)
        continue;
    end
    y_bands(:, idx) = filter(b_coeffs{idx}, a_coeffs{idx}, x);
    gain_lin = 10^(gains_dB(idx) / 20);
    y_bands(:, idx) = y_bands(:, idx) * gain_lin;
end

y_eq = sum(y_bands, 2);

% Align FIR output (constant group delay = order/2 samples).
if strcmp(filter_type, 'FIR')
    delay_samples = floor(order / 2);
    y_eq_aligned = y_eq(delay_samples+1 : end);
    x_aligned    = x(1 : end-delay_samples);
    t_aligned    = t(1 : end-delay_samples);
else
    y_eq_aligned = y_eq;
    x_aligned    = x;
    t_aligned    = t;
end

% FIX 8: Normalize output to [-1, 1] before playback/save to prevent
% clipping artefacts when gains push the signal out of range.
max_val = max(abs(y_eq_aligned));
if max_val > 1
    y_eq_aligned = y_eq_aligned / max_val;
    fprintf('\n[Note] Output was normalized to prevent clipping.\n');
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
NFFT  = length(x_aligned);
X     = abs(fftshift(fft(x_aligned)));
Y     = abs(fftshift(fft(y_eq_aligned)));
f_axis = linspace(-Fs/2, Fs/2, NFFT);

figure('Name', 'Magnitude Spectra');
plot(f_axis, 20*log10(X + eps), 'b'); hold on;   % eps avoids log10(0)
plot(f_axis, 20*log10(Y + eps), 'r');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Magnitude spectrum: original vs. equalised');
legend('Original', 'Equalised');
grid on; xlim([-Fs/2, Fs/2]);

% ----- 3. Welch Power Spectral Density (manual implementation) -----------
Nfft_welch  = 256;
overlap_frac = 0.5;
% FIX 6: Pass the number of OVERLAPPING samples (not the step size).
% noverlap = round(Nfft_welch * overlap_frac) = 128.
% Inside welch_manual, seg_shift = nfft - noverlap = 128. Correct.
noverlap_welch = round(Nfft_welch * overlap_frac);
window_welch   = hamming(Nfft_welch);   % FIX 2 applied here too

[psd_orig, f_welch] = welch_manual(x_aligned,    window_welch, Nfft_welch, noverlap_welch, Fs);
[psd_eq,   ~]       = welch_manual(y_eq_aligned, window_welch, Nfft_welch, noverlap_welch, Fs);

figure('Name', 'Welch PSD');
plot(f_welch, 10*log10(psd_orig + eps), 'b'); hold on;
plot(f_welch, 10*log10(psd_eq   + eps), 'r');
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
title('Welch PSD: original (blue) vs. equalised (red)');
legend('Original', 'Equalised');
grid on;

% ----- 4. Spectrogram (manual STFT) --------------------------------------
win_len      = 256;
overlap_stft = round(win_len * 0.75);

[S_orig, f_stft, t_stft] = stft_manual(x_aligned,    win_len, overlap_stft, Fs);
[S_eq,   ~,      ~]      = stft_manual(y_eq_aligned, win_len, overlap_stft, Fs);

figure('Name', 'Spectrograms');
subplot(2,1,1);
imagesc(t_stft, f_stft, 20*log10(abs(S_orig) + eps)); axis xy; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Original spectrogram');
ylim([0, Fs/2]);
% FIX 7: caxis is deprecated in R2022a+; use clim instead (backwards compatible via try/catch).
try, clim([-60 0]); catch, caxis([-60 0]); end

subplot(2,1,2);
imagesc(t_stft, f_stft, 20*log10(abs(S_eq) + eps)); axis xy; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Equalised spectrogram');
ylim([0, Fs/2]);
try, clim([-60 0]); catch, caxis([-60 0]); end

% ----- 5. Filter characteristics for each band --------------------------
for idx = 1:num_bands
    % Skip bands that were degenerate (above Nyquist).
    if ~band_active(idx)
        continue;
    end
    b = b_coeffs{idx};
    a = a_coeffs{idx};

    figure('Name', sprintf('Band %d characteristics', idx));

    % FIX 4 & 5: freqz() called inside a subplot creates its own figure
    % and ignores the active subplot. Capture outputs and plot manually.
    [H, f_hz] = freqz(b, a, 1024, Fs);

    subplot(2, 3, 1);
    plot(f_hz, 20*log10(abs(H) + eps));
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title(sprintf('Band %d: Magnitude Response', idx));
    grid on;

    subplot(2, 3, 2);
    plot(f_hz, unwrap(angle(H)) * (180/pi));
    xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
    title(sprintf('Band %d: Phase Response', idx));
    grid on;

    % FIX 4: impz(b, a, [], Fs) is not a valid standard call.
    % Compute the impulse response with impz(b, a) which returns samples,
    % then build a time vector manually using Fs.
    [h_imp, n_imp] = impz(b, a);
    t_imp = n_imp / Fs;   % convert sample indices to seconds

    subplot(2, 3, 3);
    stem(t_imp, h_imp, 'filled', 'MarkerSize', 3);
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Impulse Response');
    grid on;

    % Step response: cumulative sum of impulse response
    step_resp = cumsum(h_imp);
    subplot(2, 3, 4);
    plot(t_imp, step_resp);
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Step Response');
    grid on;

    % Pole-zero plot
    subplot(2, 3, 5);
    zplane(b, a);
    title('Pole-Zero Plot');
end

%% =============== PLAYBACK AND SAVE ======================================
fprintf('\nPlaying the equalised audio...\n');
sound(y_eq_aligned, Fs);
audiowrite('equalized_output.wav', y_eq_aligned, Fs);
fprintf('Equalised audio saved as "equalized_output.wav".\n');

%% =============== HELPER FUNCTIONS =======================================

function [psd, f] = welch_manual(x, window, nfft, noverlap, Fs)
% Manual Welch PSD estimate.
%   x        - input signal (column vector)
%   window   - analysis window of length nfft
%   nfft     - FFT length (= window length)
%   noverlap - number of OVERLAPPING samples between successive segments
%   Fs       - sampling frequency
    win = window(:);
    seg_shift = nfft - noverlap;   % hop size in samples
    if seg_shift <= 0
        error('Overlap too large; seg_shift must be positive.');
    end
    n_seg   = floor((length(x) - noverlap) / seg_shift);
    psd_sum = zeros(nfft, 1);
    for k = 1:n_seg
        start_idx    = (k-1)*seg_shift + 1;
        seg          = x(start_idx : start_idx + nfft - 1);
        seg_windowed = seg .* win;
        seg_fft      = fft(seg_windowed, nfft);
        psd_sum      = psd_sum + abs(seg_fft).^2;
    end
    psd = psd_sum / (n_seg * Fs * sum(win.^2));
    psd = psd(1 : floor(nfft/2)+1);
    psd(2:end-1) = 2 * psd(2:end-1);   % single-sided compensation
    f = linspace(0, Fs/2, length(psd));
end

function [S, f, t] = stft_manual(x, win_len, noverlap, Fs)
% Manual Short-Time Fourier Transform.
%   x        - input signal (column vector)
%   win_len  - window length in samples
%   noverlap - number of overlapping samples
%   Fs       - sampling frequency
%
% Returns:
%   S  - complex STFT matrix (frequency x time frames), positive freqs only
%   f  - frequency vector (Hz)
%   t  - time vector (s), centres of each analysis frame

    % FIX 2 applied here: hamming() with no extra arguments.
    win       = hamming(win_len);
    seg_shift = win_len - noverlap;
    if seg_shift <= 0
        error('Invalid overlap; seg_shift must be positive.');
    end
    nfft  = win_len;
    n_seg = floor((length(x) - noverlap) / seg_shift);
    S_full = zeros(nfft, n_seg);
    for k = 1:n_seg
        start_idx    = (k-1)*seg_shift + 1;
        seg          = x(start_idx : start_idx + win_len - 1);
        seg_windowed = seg .* win;
        S_full(:, k) = fft(seg_windowed, nfft);
    end
    % Keep only positive frequencies up to Nyquist
    keep = 1 : floor(nfft/2)+1;
    S    = S_full(keep, :);
    f    = linspace(0, Fs/2, length(keep))';
    t    = ((0:n_seg-1) * seg_shift + (win_len/2)) / Fs;
end