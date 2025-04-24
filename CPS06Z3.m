clear; close all; clc;
% Parametry
fs = 3.2e6;          % częstotliwość próbkowania
N  = 32e6;           % liczba próbek
fc = 407812.5;      % częstotliwość przesunięcia
bwSERV = 80e3;       % pasmo stacji FM
bwAUDIO = 16e3;      % pasmo audio mono

% Wczytanie surowych próbek (8-bit unsigned IQ: [I0 Q0 I1 Q1 ...])
fid = fopen('samples_100MHz_fs3200kHz.raw','rb');
s = fread(fid, 2*N, 'uint8');
fclose(fid);
s = int16(s) - 127;

% Rozdział na I i Q oraz złożenie sygnału zespolonego
I = double(s(1:2:end));
Q = double(s(2:2:end));
wideband = I + 1j*Q;
t = (0:N-1)'/fs;

% Shift do zera stacji
wideband_shifted = wideband .* exp(-1j*2*pi*fc*t);

% Projektowanie i aplikacja filtru dolnoprzepustowego (80 kHz)
[b_lp,a_lp] = butter(4, bwSERV/(fs/2), 'low');
wb_filt = filter(b_lp, a_lp, wideband_shifted);

% Decymacja do pasma ±bwSERV
dec_factor = floor(fs/(2*bwSERV));
x = wb_filt(1:dec_factor:end);

% FM-demodulacja przez różnicę fazy
dx = x(2:end) .* conj(x(1:end-1));
y = angle(dx);

% Filtr dolnoprzepustowy audio (16 kHz przy fs=160kHz)
fs1 = 160e3;
[b2,a2] = butter(4, bwAUDIO/(fs1/2), 'low');
y_filt = filter(b2, a2, y);

% Decymacja do 32 kHz
dec2 = floor(fs1/(2*bwAUDIO));
ym = y_filt(1:dec2:end);

% Filtr de-emfazy (τ = 75 µs) przy fs_audio = 32 kHz
fs_audio = 32e3;
tau = 75e-6;
dt = 1/fs_audio;
alpha = dt/(tau+dt);
b_de = alpha;
a_de = [1, alpha-1];
ym = filter(b_de, a_de, ym);

% Normalizacja sygnału audio
ym = ym - mean(ym);
ym = ym / (1.001 * max(abs(ym)));
soundsc( ym, fs_audio);


% --- PSD (Pwelch) ---
figure; 
[pxx,f] = pwelch(wideband, 2048, [], [], fs);
semilogy(f,pxx); grid on;
title('|X(f)| wideband'); xlabel('Hz'); ylabel('PSD');

figure; 
[pxx,f] = pwelch(wideband_shifted, 2048, [], [], fs);
semilogy(f,pxx); grid on;
title('|X(f)| shifted'); xlabel('Hz'); ylabel('PSD');

figure; 
[pxx,f] = pwelch(wb_filt, 2048, [], [], fs);
semilogy(f,pxx); grid on;
title('|X(f)| filtered'); xlabel('Hz'); ylabel('PSD');

figure; 
[pxx,f] = pwelch(x, 2048, [], [], fs);
semilogy(f,pxx); grid on;
title('|X(f)| decimated'); xlabel('Hz'); ylabel('PSD');

figure; 
[pxx,f] = pwelch(y, 2048, [], [], fs1);
semilogy(f,pxx); grid on;
title('|X(f)| demodulated'); xlabel('Hz'); ylabel('PSD');

figure; 
[pxx,f] = pwelch(ym, 2048, [], [], fs_audio);
semilogy(f,pxx); grid on;
title('|X(f)| ym decimated down sampled'); xlabel('Hz'); ylabel('PSD');


% wykres przebiegu czasowego mono-audio
fs_audio   = 32e3;            % już zdefiniowane w Twoim skrypcie
N_audio    = length(ym);      % liczba próbek po wszystkich operacjach
t_audio    = (0:N_audio-1)/fs_audio;   % wektor czasu w [s]
figure;
plot(t_audio, ym);
xlabel('t  [s]');
ylabel('amplitude');
title('decoded / down-sampled mono audio');
grid on;
xlim([0 10]);
ylim([-2 2]);    



% % --- Spektrogramy ---
% figure;
% spectrogram(wideband_shifted, 1024, 512, 1024, fs, 'yaxis');
% title('Spectrogram shifted');
% 
% figure;
% spectrogram(wb_filt, 1024, 512, 1024, fs, 'yaxis');
% title('Spectrogram filtered');
% 
% figure;
% spectrogram(y_filt, 1024, 512, 1024, fs1, 'yaxis');
% title('Spectrogram y\_filt');
% 
% figure;
% spectrogram(ym, 1024, 512, 1024, fs_audio, 'yaxis');
% title('Spectrogram audio');


%% --- Spektrogramy ---  (ustawienia okno = 1024, overlap = 3/4)
wlen  = 1024;
nov   = wlen/2;
nfft  = wlen;

% 1. oryginalny strumień I/Q (IF ≈ 100 kHz)
figure;
spectrogram(wideband, wlen, nov, nfft, fs, 'yaxis');
title('Spectrogram ①  wideband  (Fs = 3.2 MHz)');

% 2. po przesunięciu do 0 Hz
figure;
spectrogram(wideband_shifted, wlen, nov, nfft, fs, 'yaxis');
title('Spectrogram ②  shifted to baseband');

% 3. po filtrze LP 80 kHz
figure;
spectrogram(wb_filt, wlen, nov, nfft, fs, 'yaxis');
title('Spectrogram ③  LP 80 kHz');

% 4. po decymacji ×20  →  Fs = 160 kHz
figure;
spectrogram(x, wlen, nov, nfft, fs1, 'yaxis');
title('Spectrogram ④  x  (Fs = 160 kHz)');

% 5. „sygnał hybrydowy” tuż po demodulacji FM
figure;
spectrogram(y, wlen, nov, nfft, fs1, 'yaxis');
title('Spectrogram ⑤  y  – demodulated hybrid (mono + 19 kHz pilot + stereo + RDS)');

% 6. zdekodowany mono-audio po LP 16 kHz, decymacji ×5 i de-emfazie
figure;
spectrogram(ym, wlen, nov, nfft, fs_audio, 'yaxis');
title('Spectrogram ⑥  ym  (mono audio, Fs = 32 kHz)');




% --- Charakterystyki częstotliwościowe filtrów de- and pre-emfazy ---
f_cutoff = 2100;  % Hz
order = 1;

% de-emphasis (mały LP przy 2.1 kHz)
[b_de2,a_de2] = butter(order, f_cutoff/(fs_audio/2), 'low');
% pre-emphasis (HP przy 2.1 kHz)
[b_pre,a_pre] = butter(order, f_cutoff/(fs_audio/2), 'high');

[H_de,w] = freqz(b_de2, a_de2, 8000, fs_audio);
[H_pr,~] = freqz(b_pre, a_pre, 8000, fs_audio);

figure; 
semilogx(w,20*log10(abs(H_de)+eps)); hold on;
semilogx(w,20*log10(abs(H_pr)+eps));
xline(f_cutoff,'--g','2.1 kHz');
grid on; legend('de-emphasis','pre-emphasis');
title('Charakterystyki amplitudowo-częstotliwościowe');
xlabel('Hz'); ylabel('dB');
