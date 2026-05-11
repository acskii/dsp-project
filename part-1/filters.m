%%
fs = 360;           % Sampling frequency
fc = 0.5;           % Cutoff frequency
order = 4;          % Filter order

[b, a] = butter(order, fc/(fs/2), 'high');

figure;
freqz(b, a, 1024, fs); 
title('IIR Butterworth HPF: Magnitude & Phase');

figure;
subplot(2,2,1); 
zplane(b, a); 
title('Pole-Zero');

subplot(2,2,2); 
impz(b, a, 100, fs); 
title('Impulse Response');

subplot(2,2,3); 
stepz(b, a, 100, fs); 
title('Step Response');

%%
fs = 360;           % Sampling frequency
f_notch = 50;       % Frequency of the notch
bw = 2;             % Bandwidth of the notch
[b, a] = iirnotch(f_notch/(fs/2), bw/(fs/2));

% Analysis
figure;
freqz(b, a, 1024, fs); 
title('Notch Filter: Magnitude & Phase');

figure;
subplot(2,2,1); 
zplane(b, a); 
title('Pole-Zero');

subplot(2,2,2); 
impz(b, a, 100, fs); 
title('Impulse Response');

subplot(2,2,3); 
stepz(b, a, 100, fs); 
title('Step Response');

%%
fs = 360;           % Sampling frequency
fc = 100;           % Cutoff for Muscle Noise
num_taps = 51;      % 50 + 1
b = fir1(num_taps - 1, fc/(fs/2), 'low', hamming(num_taps));

% Analysis
figure;
freqz(b, 1, 1024, fs); 
title('FIR Window-based Filter: Magnitude & Phase');

figure;
subplot(2,2,1); 
zplane(b, 1); 
title('Pole-Zero');

subplot(2,2,2); 
stem(b); 
title('Impulse Response');

subplot(2,2,3); 
stepz(b, 1, 100, fs); 
title('Step Response');
