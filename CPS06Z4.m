clear; close all; clc;
[y1, fs1] = audioread('motor.wav');
[y2, fs2] = audioread('bird.wav');

% Dopasowanie 
min_len = min(length(y1), length(y2));
y1 = y1(1:min_len);
y2 = y2(1:min_len);
y_total = y1 + y2;

% Parametry FFT
N_total = length(y_total);
% FFT sumy
Y_total = fft(y_total);
f = (0:N_total-1)*(fs1/N_total);

figure;
plot(f(1:N_total/2), abs(Y_total(1:N_total/2)));
title('Widmo FFT suma'); xlabel('Czesntotliwosc [Hz]'); ylabel('Amplituda'); grid on;

% FFT silnika
N1 = length(y1);
Y1 = fft(y1);
f1 = (0:N1-1)*(fs1/N1);

figure;
plot(f1(1:N1/2), abs(Y1(1:N1/2)));
title('Widmo FFT - silnik'); xlabel('Czesntotliwosc [Hz]'); ylabel('Amplituda'); grid on;

% FFT ptaka
N2 = length(y2);
Y2 = fft(y2);
f2 = (0:N2-1)*(fs2/N2);

figure;
plot(f2(1:N2/2), abs(Y2(1:N2/2)));
title('Widmo FFT - ptak'); xlabel('Czesntotliwosc [Hz]'); ylabel('Amplituda'); grid on;

% Spektrogramy
window = 1024;
noverlap = 512;
nfft = 1024;

figure;
spectrogram(y_total, window, noverlap, nfft, fs1, 'yaxis');
title('Spektrogram - suma');

figure;
spectrogram(y1, window, noverlap, nfft, fs1, 'yaxis');
title('Spektrogram - silnik');

figure;
spectrogram(y2, window, noverlap, nfft, fs2, 'yaxis');
title('Spektrogram - ptak');

% Filtr dolnoprzepustowy Butterworth
cutoff = 1500; % Hz
order = 4;
[b, a] = butter(order, cutoff/(fs1/2), 'low');
y_filtered = filter(b, a, y_total);

% FFT po filtrze
Nf = length(y_filtered);
Yf = fft(y_filtered);
ff = (0:Nf-1)*(fs1/Nf);

figure;
plot(ff(1:Nf/2), abs(Yf(1:Nf/2)));
title('Widmo FFT - po filtrze'); xlabel('Czesntotliwosc [Hz]'); ylabel('Amplituda'); grid on;

% Spektrogram po filtrze
figure;
spectrogram(y_filtered, window, noverlap, nfft, fs1, 'yaxis');
title('Spektrogram - po filtrze');

% Bieguny i zera filtru
figure;
zplane(b, a);
title('Bieguny i zera filtru Butterworth');
