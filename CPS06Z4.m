%% Zadanie 4 – filtracja dźwięków rzeczywistych (wersja: motor + bird, poprawka zakresu filtru) 
% (c) Twoje Imię 2024 – wolno korzystać i modyfikować
% -------------------------------------------------------------------------
%  W tej wersji nie korzystamy z pliku „speech.wav”.  Miks składa się z:
%  • motor.wav – niskoczęstotliwościowy warkot silnika
%  • bird.wav  – wysokoczęstotliwościowy śpiew ptaka
%  Celem jest odseparowanie śpiewu ptaka od tła silnika.
% -------------------------------------------------------------------------

clear; close all; clc;

%% 1. Wczytanie sygnałów ---------------------------------------------------
[motor, fsMotor] = audioread("10_14B.wav");   % burza
[bird , fsBird ] = audioread("grilli.wav");    % świerszcze (dawniej ptak)

%% 2. Wyrównanie częstotliwości próbkowania --------------------------------
fs   = max([fsMotor, fsBird]);                % użyj najwyższej dostępnej
motor = resample(motor, fs, fsMotor);
bird  = resample(bird , fs, fsBird );

%% 3. Ucięcie do wspólnej długości i stworzenie sumy -----------------------
L     = min([numel(motor), numel(bird)]);
motor = motor(1:L);
bird  = bird (1:L);

mix   = motor + bird;
mix   = mix ./ max(abs(mix));                 % normalizacja – zapobiega przesterowi

%% 4. Widma FFT pojedynczych sygnałów i sumy --------------------------------
figure('Name','Widma FFT','Position',[100 100 1400 600])
subplot(2,2,1);  plotFFT(motor, fs, 'Motor');
subplot(2,2,2);  plotFFT(bird , fs, 'Bird');
subplot(2,2,3);  plotFFT(mix  , fs, 'Motor + Bird');

%% 5. Spektrogramy (STFT) ---------------------------------------------------
figure('Name','Spektrogramy');
subplot(2,2,1);  plotSpec(motor, fs, 'Motor');
subplot(2,2,2);  plotSpec(bird , fs, 'Bird');
subplot(2,2,3);  plotSpec(mix  , fs, 'Sum (mix)');

%% 6. Projekt rekursywnego filtru IIR – zostawiamy tylko śpiew ptaka --------
%  Zakładamy, że śpiew ptaka występuje głównie > 3 kHz, a warkot silnika
%  poniżej 1 kHz.  Poniższy kod sam dopasowuje górną granicę pasma do
%  aktualnego fs, żeby uniknąć błędu „frequencies must be over the (0,1)
%  interval” przy zbyt niskiej częstotliwości próbkowania.
order          = 6;
lowCutDesired  = 150;        % Hz
highCutDesired = 2000;       % Hz (domyślnie)
nyq            = fs/2;        % częst. Nyquista

% Skoryguj górną granicę, jeśli wykracza poza Nyquista
highCut = min(highCutDesired, 0.9*nyq);   % 0.9*nyq zostawia zapas
lowCut  = lowCutDesired;

% Jeśli pasmo „zwraca się” – przejdź na filtr górnoprzepustowy
if lowCut >= highCut
    warning("fs=%.0f Hz jest zbyt małe na pasmo 3–12 kHz – używam filtru HPF > %.0f Hz", fs, lowCut);
    [b,a] = butter(order, lowCut/nyq, 'high');
    filtType = sprintf('High‑pass > %d Hz', lowCut);
else
    [b,a] = butter(order, [lowCut highCut]/nyq, 'bandpass');
    filtType = sprintf('Band‑pass %d–%d Hz', lowCut, round(highCut));
end

% Charakterystyka amplitudowa w dB ----------------------------------------
figure('Name','Charakterystyka filtru');
[H,w] = freqz(b, a, 4096, fs);
plot(w, 20*log10(abs(H))), grid on
xlabel('Częstotliwość [Hz]'); ylabel('|H(f)| [dB]')
ylim([-80 5]); title(['Butterworth ' filtType])

% Zera i bieguny -----------------------------------------------------------
figure('Name','Z-P filtru');  zplane(b, a);  title(['Z‑plane ' filtType])

%% 7. Filtracja sygnału sumy i odsłuch --------------------------------------
filtered = filter(b, a, mix);
filtered = filtered ./ max(abs(filtered));     % normalizacja
soundsc(filtered, fs)                         % odsłuch

% Zapis do pliku
outfile = "bird_isolated.wav";
audiowrite(outfile, filtered, fs);

%% 8. Widma FFT + STFT po filtrze ------------------------------------------
figure('Name','Widmo sygnału po filtrze');
plotFFT(filtered, fs, sprintf('Filtered (%s)', filtType));

figure('Name','Spektrogram po filtrze');
plotSpec(filtered, fs, sprintf('Filtered (%s)', filtType));

%% 9. Porównawcze zestawienie (bird vs. po filtrze) -------------------------
figure('Name','Porównanie widm – oryginalny bird vs. wynik filtru');
[Bw, faxis] = getFFT(bird    , fs);
[Fw, ~   ]  = getFFT(filtered, fs);
plot(faxis, Bw, 'DisplayName','Oryginalny bird'); hold on
plot(faxis, Fw, 'DisplayName','Po filtrze');
legend show; grid on; xlim([0 nyq])
xlabel('Częstotliwość [Hz]'); ylabel('|X(f)|');
title('Porównanie FFT');

disp("✓ Gotowe! Izolowany śpiew ptaka zapisano w pliku " + outfile);

%% -------------------------------------------------------------------------
%  Lokalne funkcje pomocnicze ----------------------------------------------
function plotFFT(x, fs, ttl)
    [X,f] = getFFT(x, fs);
    plot(f, X); grid on
    title(['FFT – ' ttl]); xlabel('Częstotliwość [Hz]'); ylabel('|X(f)|');
    xlim([0 fs/2]);
end

function [X,f] = getFFT(x, fs)
    N  = numel(x);
    win = hann(N);
    Xc = fft(x .* win);
    X  = abs(Xc(1:floor(N/2))) / N * 2;
    f  = (0:numel(X)-1)' * fs/N;
end

function plotSpec(x, fs, ttl)
    window   = 1024;
    noverlap = 512;
    nfft     = 1024;
    spectrogram(x, window, noverlap, nfft, fs, 'yaxis');
    title(['Spektrogram – ' ttl]);
end
