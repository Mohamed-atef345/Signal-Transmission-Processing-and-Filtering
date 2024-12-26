%--------------------------------------------------------------------------
% Signal Transmission, Processing, and Filtering (communication system simulation)
%--------------------------------------------------------------------------
% This script reads an audio file, processes it through various channels, 
% adds noise, filters the signal, and plays the processed audio back. 
% Visualizations are provided in both time and frequency domains.

%------------------------------ 1. Transmission ---------------------------
% Read the sound file and play it
[x, fs] = audioread('audio_file.mp3'); % Load audio file
disp('Now playing your audio file ...');
value = 1;
while (value == 1)
    sound(x, fs); % Play audio
    value = input('Press Zero to stop the audio file: '); % Prompt to stop
    if (value == 0)  
        clear sound; % Stop audio playback
    end
end

% Get the number of audio channels
channels_num = size(x, 2);

% Plot the signal in the time domain
t = linspace(0, length(x) / fs, length(x));
figure; subplot(2, 2, [1, 2]);
plot(t, x);
title('Signal in Time Domain'); xlabel('Time (s)'); ylabel('Amplitude');

% Plot the signal in the frequency domain
X = fftshift(fft(x)); % Perform FFT
fvec = linspace(-fs / 2, fs / 2, length(X)); % Frequency vector
Xmag = abs(X); % Magnitude spectrum
Xphase = angle(X); % Phase spectrum
subplot(2, 2, 3);
plot(fvec, Xmag);
title('Magnitude in Frequency Domain'); xlabel('Frequency (Hz)'); ylabel('Magnitude');
subplot(2, 2, 4);
plot(fvec, Xphase);
title('Phase in Frequency Domain'); xlabel('Frequency (Hz)'); ylabel('Phase (rad)');

%------------------------------- 2. Channel -------------------------------
% Define the channels
h1 = [1 zeros(1, length(x) - 1)]; % Channel 1: Direct pass
h2 = exp(-2 * pi * 5000 * t);     % Channel 2: Exponential decay at 5000 Hz
h3 = exp(-2 * pi * 1000 * t);     % Channel 3: Exponential decay at 1000 Hz
h4 = [2 zeros(1, (1 - 0) * fs - 2) 0.5 zeros(1, length(x) - (1 - 0) * fs)]; % Channel 4: Mixed response

% Preallocate arrays for channel outputs
x_channel_1 = zeros(length(x) * 2 - 1, channels_num);
x_channel_2 = zeros(length(x) * 2 - 1, channels_num);
x_channel_3 = zeros(length(x) * 2 - 1, channels_num);
x_channel_4 = zeros(length(x) * 2 - 1, channels_num);

% Convolve signal with each channel
for i = 1:channels_num
    x_channel_1(:, i) = conv(x(:, i), h1'); % Convolution for Channel 1
    x_channel_2(:, i) = conv(x(:, i), h2'); % Convolution for Channel 2
    x_channel_3(:, i) = conv(x(:, i), h3'); % Convolution for Channel 3
    x_channel_4(:, i) = conv(x(:, i), h4'); % Convolution for Channel 4
end

% Plot signals after passing through channels
t_conv = linspace(0, (length(x) / fs) * 2, length(x) * 2 - 1);
figure;
subplot(2, 2, 1); plot(t, x_channel_1(1:length(x), :)); title('Signal after Channel 1');
subplot(2, 2, 2); plot(t, x_channel_2(1:length(x), :)); title('Signal after Channel 2');
subplot(2, 2, 3); plot(t, x_channel_3(1:length(x), :)); title('Signal after Channel 3');
subplot(2, 2, 4); plot(t, x_channel_4(1:length(x), :)); title('Signal after Channel 4');

%-------------------------------- 3. Noise --------------------------------
% Generate noise
sigma = input('Enter the value of sigma for noise generation (Z(t) = sigma*randn(1,length(x)):\n');
Z = sigma * randn(1, length(x_channel_1)); % Random noise

% Add noise to each channel's output
x_noise_channel_1 = x_channel_1 + Z';
x_noise_channel_2 = x_channel_2 + Z';
x_noise_channel_3 = x_channel_3 + Z';
x_noise_channel_4 = x_channel_4 + Z';

% Plot noisy signals
figure;
subplot(2, 2, 1); plot(t, x_noise_channel_1(1:length(x), :)); title('Signal on Channel 1 with Noise');
subplot(2, 2, 2); plot(t, x_noise_channel_2(1:length(x), :)); title('Signal on Channel 2 with Noise');
subplot(2, 2, 3); plot(t, x_noise_channel_3(1:length(x), :)); title('Signal on Channel 3 with Noise');
subplot(2, 2, 4); plot(t, x_noise_channel_4(1:length(x), :)); title('Signal on Channel 4 with Noise');

%-------------------------------- 4. Receiver -----------------------------
% Transform noisy signals to frequency domain
X_filtered_channel_1 = fftshift(fft(x_noise_channel_1(1:length(x), :)));
X_filtered_channel_2 = fftshift(fft(x_noise_channel_2(1:length(x), :)));
X_filtered_channel_3 = fftshift(fft(x_noise_channel_3(1:length(x), :)));
X_filtered_channel_4 = fftshift(fft(x_noise_channel_4(1:length(x), :)));

% Define a low-pass filter
low_critical_point = -3400; % Lower frequency cutoff
high_critical_point = 3400; % Upper frequency cutoff
filter = zeros(size(X_filtered_channel_1)); % Initialize filter
filter(fvec > -3400 & fvec < 3400, :) = 1; % Apply frequency range

% Apply filter to each channel's spectrum
X_filtered_channel_1 = X_filtered_channel_1 .* filter;
X_filtered_channel_2 = X_filtered_channel_2 .* filter;
X_filtered_channel_3 = X_filtered_channel_3 .* filter;
X_filtered_channel_4 = X_filtered_channel_4 .* filter;

% Convert filtered signals back to time domain
x_filtered_channel_1 = real(ifft(ifftshift(X_filtered_channel_1)));
x_filtered_channel_2 = real(ifft(ifftshift(X_filtered_channel_2)));
x_filtered_channel_3 = real(ifft(ifftshift(X_filtered_channel_3)));
x_filtered_channel_4 = real(ifft(ifftshift(X_filtered_channel_4)));

% Plot filtered signals in time domain
figure;
subplot(2, 2, 1); plot(t, x_filtered_channel_1); title('Filtered Signal from Channel 1 in Time Domain');
subplot(2, 2, 2); plot(t, x_filtered_channel_2); title('Filtered Signal from Channel 2 in Time Domain');
subplot(2, 2, 3); plot(t, x_filtered_channel_3); title('Filtered Signal from Channel 3 in Time Domain');
subplot(2, 2, 4); plot(t, x_filtered_channel_4); title('Filtered Signal from Channel 4 in Time Domain');

% Plot filtered signals in frequency domain
figure;
subplot(2, 2, 1); plot(fvec, abs(X_filtered_channel_1)); title('Filtered Signal from Channel 1 in Frequency Domain');
subplot(2, 2, 2); plot(fvec, abs(X_filtered_channel_2)); title('Filtered Signal from Channel 2 in Frequency Domain');
subplot(2, 2, 3); plot(fvec, abs(X_filtered_channel_3)); title('Filtered Signal from Channel 3 in Frequency Domain');
subplot(2, 2, 4); plot(fvec, abs(X_filtered_channel_4)); title('Filtered Signal from Channel 4 in Frequency Domain');

% Play filtered signals
disp('Now playing filtered signal from Channel 1 ...');
sound(x_filtered_channel_1, fs);
pause(length(x) / fs);
disp('Now playing filtered signal from Channel 2 ...');
sound(x_filtered_channel_2, fs);
pause(length(x) / fs);
disp('Now playing filtered signal from Channel 3 ...');
sound(x_filtered_channel_3, fs);
pause(length(x) / fs);
disp('Now playing filtered signal from Channel 4 ...');
sound(x_filtered_channel_4, fs);
pause(length(x) / fs);
