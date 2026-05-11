%%
fs = 360;            % Sampling frequency
num_samples = 5 * fs;

% logic to read binary Format 212
records = {'100', '106'};
data_signals = cell(1, 2);

for i = 1:2
    fname = [records{i} '.dat'];
    fid = fopen(fname, 'r');
    % MIT-BIH Format 212: 3 bytes represent 2 samples (1 from each channel)
    A = fread(fid, [3, num_samples], 'uint8')';
    fclose(fid);
    
    % Bit-shifting logic to decode Format 212
    M = bitshift(bitand(A(:,2), 15), 8) + A(:,1);
    % Convert to mV (approximate gain for these records is 200)
    data_signals{i} = (double(M) - 1024) / 200; 
end

ecg100 = data_signals{1};
ecg106 = data_signals{2};
t = (0:length(ecg100)-1) / fs;

% High Pass Filter
fc1 = 0.5;           % Cutoff frequency
order = 4;           % Filter order
[b_hp, a_hp] = butter(order, fc1/(fs/2), 'high');

% Notch Filter
f_notch = 50;       % Frequency of the notch
bw = 2;             % Bandwidth of the notch
[b_n, a_n] = iirnotch(f_notch/(fs/2), bw/(fs/2));

% Low Pass Filter
fc2 = 100;           % Cutoff for Muscle Noise
num_taps = 51;       % 50 + 1
b_lp = fir1(num_taps - 1, fc2/(fs/2), 'low', hamming(num_taps));

% Filtering signal
% Record 100
ecg100_1    = filtfilt(b_hp, a_hp, ecg100);
ecg100_2 = filtfilt(b_n, a_n, ecg100_1);
ecg100_final  = filtfilt(b_lp, 1, ecg100_2);

% Record 106
ecg106_1    = filtfilt(b_hp, a_hp, ecg106);
ecg106_2 = filtfilt(b_n, a_n, ecg106_1);
ecg106_final  = filtfilt(b_lp, 1, ecg106_2);

% Comparison
figure;
subplot(2,2,1); 
plot(t, ecg100, 'r');
title('Record 100 - Raw'); 
xlabel('Time (s)'); 
ylabel('mV'); 
grid on;

subplot(2,2,2); 
plot(t, ecg100_final, 'b');
title('Record 100 - Filtered'); 
xlabel('Time (s)'); 
ylabel('mV'); 
grid on;

subplot(2,2,3); 
plot(t, ecg106, 'r');
title('Record 106 - Raw'); 
xlabel('Time (s)'); 
ylabel('mV'); 
grid on;

subplot(2,2,4); 
plot(t, ecg106_final, 'b');
title('Record 106 - Filtered'); 
xlabel('Time (s)'); 
ylabel('mV'); 
grid on;

% Power spectral density
figure;
subplot(2,1,1);
pwelch(ecg100, [], [], [], fs); 
hold on;
pwelch(ecg100_final, [], [], [], fs);
title('PSD Comparison - Record 100');
legend('Raw', 'Filtered');

subplot(2,1,2);
pwelch(ecg106, [], [], [], fs); 
hold on;
pwelch(ecg106_final, [], [], [], fs);
title('PSD Comparison - Record 106');
legend('Raw', 'Filtered');

% Spectrogram
figure;
subplot(2,1,1);
spectrogram(ecg100_final, hamming(128), 120, 128, fs, 'yaxis');
title('Spectrogram of Filtered Record 100');

subplot(2,1,2);
spectrogram(ecg106_final, hamming(128), 120, 128, fs, 'yaxis');
title('Spectrogram of Filtered Record 106');

% SNR Improvement
noise100 = ecg100 - ecg100_final;
noise106 = ecg106 - ecg106_final;
snr100 = 10 * log10(var(ecg100_final) / var(noise100));
snr106 = 10 * log10(var(ecg106_final) / var(noise106));
fprintf('SNR Improvement Summary: \n');
fprintf('Record 100 : %.2f dB\n', snr100);
fprintf('Record 106 : %.2f dB\n', snr106);