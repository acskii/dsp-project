function DSP_Final()
    audioName = input('Enter the audio filename: ', 's');
    [signal, fs] = loadAudio(audioName);
    
    bands = chooseMode();
    
    gains = inputGains(bands);
end

function [signal, fs] = loadAudio(fileName)
    % Input: 
    %   fileName: String of file name within code directory
    %
    % Output:
    %   signal: Mono audio data normalised between -1 and 1
    %   fs: Sampling frequency
    
    try 
        [raw, fs] = audioread(fileName);
        
        if size(raw, 2) > 1
            fprintf('Stereo signal detected. Converting to mono...\n');
            signal = mean(raw, 2);
        else
            signal = raw;
        end
        
        % Normalisation
        signal = signal / max(abs(signal));
        
        fprintf('Successfully loaded: %s\n', fileName);
        fprintf('Sampling Rate: %d Hz | Duration: %.2f seconds\n', fs, length(signal)/fs);
    catch err
       error('Error loading audio file: %s.', err.message); 
    end
end

function bands = chooseMode()
    % Output: 
    %   bands: [f0, f1, f2, ... fn]
    
    choice = 0;
    
    while choice < 1 || choice > 2
        fprintf('\n--- Equalizer Mode Selection ---\n');
        fprintf('1. Preset Mode\n');
        fprintf('2. Custom Mode\n');
        choice = input('Select mode (1 or 2): ');
    end
    
    if choice == 1
        % Preset Mode
        bands = [0, 100, 300, 800, 2000, 5000, 10000, 20000];
    else
        % Custome Mode
        n = input('Enter number of bands (between 5 and 10): ');
        while n < 5 || n > 10
            n = input('Invalid. Please enter a number between 5 and 10: ');
        end
        
        bands = zeros(1, n);
        bands(1) = 0;
        bands(end) = 20000;
        
        for i = 2 : n - 1
            msg = sprintf('Enter end frequency for Band %d in Hz (Must be > %d): ', i - 1, bands(i - 1));
            bands(i) = input(msg);
            
            % Validation: Ensure frequencies are increasing
            while bands(i) <= bands(i - 1)
                bands(i) = input('Error: Frequency must be higher than previous band. Try again: ');
            end
        end
    end
end

function gains = inputGains(bands)
    % Prompts user for dB gain for each frequency band
    % Input: 
    %   bands [output of chooseMode()]
    % Output: 
    %   gains: [g0, g1, g2, ... gn]

    numBands = length(bands) - 1;
    gains = zeros(1, numBands);
    
    fprintf('\n--- Gain Configuration (in dB) ---\n');
    for i = 1:numBands
        prompt = sprintf('Band %d (%g Hz - %g Hz): ', i, bands(i), bands(i+1));
        gains(i) = input(prompt);
    end
end

function doBandAnalysis(b, a, fs, index)
    % Generates and saves required plots for a specific filter band
    % Does ONE BAND at a time
    % Use only after filter design and after applying gains
    %
    % Inputs:
    %   b: Cofficient of numerator of filter 
    %   a: Cofficient of denomenator of filter 
    %   fs: Sampling rate
    %   index: Band index
    % Output:
    %   Saves plots as images in filter_analysis/band{index}/
    
    % create folders
    mainDir = 'filter_analysis';
    subDir = fullfile(mainDir, sprintf('band%d', index));
    
    if ~exist(subDir, 'dir')
        mkdir(subDir); 
    end
    fig = figure('Visible', 'off');
    
    % Magnitude & Phase Response
    [H, F] = freqz(b, a, fs);
    magH = abs(H);
    phaseH = angle(H)*180/pi;
     
    plot(F, magH); 
    title(sprintf('Band %d: Magnitude Response', index));
    ylabel('Gain (dB)'); 
    grid on;
    saveas(fig, fullfile(subDir, 'magnitude.png'));
    clf(fig);
    
    plot(F, phaseH); 
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)'); 
    grid on;
    saveas(fig, fullfile(subDir, 'phase.png'));
    clf(fig);
    
    % Impulse Response
    impz(b, a, 100, fs);
    title('Impulse Response');
    saveas(fig, fullfile(subDir, 'impulse.png'));
    clf(fig);
    
    % Step Response
    stepz(b, a, 100, fs); 
    title('Step Response'); 
    xlabel('Time (s)');
    saveas(fig, fullfile(subDir, 'step.png'));
    clf(fig);
    
    % Pole-Zero
    zplane(b, a); 
    title(sprintf('Pole-Zero Plot (Order: %d)', order));
    saveas(fig, fullfile(subDir, 'zero.png'));
    clf(fig);
    
    close(fig);
end