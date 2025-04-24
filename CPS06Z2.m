%% Zadanie 2 – Dekodowanie DTMF
% -------------------------------------------------------------------------
% Skrypt dekoduje sekwencję DTMF zapisaną w plikach s0.wav…s9.wav (lub s.wav)
% oraz ilustruje przebieg działań: spektrogram oryginalny, spektrogram po
% filtracji wąskopasmowej (BP z ćw. 1) + wykres sygnał vs sygnał po filtracji.
% -------------------------------------------------------------------------
clear all; close all; clc;
%% 1. Wybór plików z danymi
fnamee = 's.wav';
[sxx, fss] = audioread(fnamee);

fname = 's4.wav'; 
[sx, fs] = audioread(fname);

%% 2. Parametry spektrogramu
win     = 4096;                         % długość okna
noverlap= win-512;                      % krok 512 samp (≈32 ms)
ff      = 0:5:2000;                     % wektor freq (Hz)

figure('Name','Spektrogram – s');
[sPP,FF,TT] = spectrogram(sxx, win, noverlap, ff, fss);
imagesc(TT, FF, 20*log10(abs(sPP)+eps)); axis xy; colormap jet;
colorbar; xlabel('t [s]'); ylabel('f [Hz]'); title('Spektrogram s [1,2,3,4,5,6,7,8,9,*,0,#]');

figure('Name','Spektrogram – s4');
[sP,F,T] = spectrogram(sx, win, noverlap, ff, fs);
imagesc(T, F, 20*log10(abs(sP)+eps)); axis xy; colormap jet;
colorbar; xlabel('t [s]'); ylabel('f [Hz]'); title('Spektrogram s4 [3 8 2 9 2]');


%% 4. Filtr BP z ćwiczenia 1 (1189–1229 Hz)
fc1 = 1189; fc2 = 1229;
[z, p, k] = butter(4, [fc1 fc2]*2*pi, 'bandpass', 's');  % analog
[ba, aa]  = zp2tf(z, p, k);
[bd, ad]  = bilinear(ba, aa, fs);

% filtracja (filtfilt – zero‑phase, więc brak przesunięcia)
sy = filtfilt(bd, ad, sx);

%% 5. Spektrogram – po filtracji
figure('Name','Spektrogram – po filtracji');
[sPf,Ff,Tf] = spectrogram(sy, win, noverlap, ff, fs, 'yaxis');
imagesc(Tf, Ff, 20*log10(abs(sPf)+eps)); axis xy; colormap jet;
colorbar; xlabel('t [s]'); ylabel('f [Hz]'); title('Po filtracji BP 1189–1229 Hz');

%% 6. Porównanie przebiegów w czasie
figure('Name','Czas – oryginał vs po BP');
plot((0:numel(sx)-1)/fs, sx); hold on;
plot((0:numel(sy)-1)/fs, sy, 'LineWidth',1.1);
legend({'oryginał','po BP'}); xlabel('t [s]'); grid on;
title(sprintf('Porównanie – %s', fname));

% Kompensacja opóźnienia (gdybyśmy użyli filter, a nie filtfilt)
% Jeśli zamiast filtfilt() używasz filter(), przesuń sygnał o d=(max(M,N)-1)/2.
% M = length(bd); N = length(ad);

%% 7. Opcjonalne (+0,25) – Transformata DtFT (algorytm Goertzla)
% -------------------------------------------------------------------------
% Dla każdego 1-sekundowego bloku obliczamy wartość DtFT w 7 interesujących
% częstotliwościach DTMF i wyświetlamy ich energie.  Bez (jeszcze) żadnej
% decyzji, filtrów ani mapowania na cyfry – to należy do kolejnych zadań.
% -------------------------------------------------------------------------
blockDur = 1.0;                         % długość bloku [s]
N        = round(blockDur*fs);          % próbki/blok
freqs    = [697 770 852 941 1209 1336 1477];   % 7 tonów DTMF
nBlocks  = floor(length(sx)/N);

energies = zeros(nBlocks, numel(freqs));

for b = 1:nBlocks
    seg = sx((b-1)*N + (1:N));         % próbki bieżącego bloku
    for k = 1:numel(freqs)
        bin        = round(freqs(k)*N/fs);  % pozycja w widmie
        energies(b,k) = abs(goertzel(seg, bin));
    end
end

% Podgląd energii (bar-plot) dla pierwszego bloku
figure('Name','Goertzel – blok #1');
stem(freqs, energies(1,:), 'filled');
xlabel('f [Hz]'); ylabel('|X(e^{jω})|'); grid on;
title('Energie DtFT (Goertzel) – blok 1');

% Dodatkowo: „mapa ciepła” blok × częstotliwość (łatwo dostrzec pary tonów)
figure('Name','Goertzel – heat-map');
imagesc(freqs, 1:nBlocks, 20*log10(energies+eps));
xlabel('f [Hz]'); ylabel('blok'); axis xy; colormap jet; colorbar;
title('Energie Goertzla w kolejnych blokach');


%% 8. Opcjonalne (+0,25) – IIR
dtmf_freqs_low = [697, 770, 852, 941];  % Hz, grupa niska
dtmf_freqs_high = [1209, 1336, 1477];   % Hz, grupa wysoka
window_duration = 0.5;  % sekundy
step_duration = 1.0;    % sekundy
threshold = 0.5;        % próg detekcji (dostosuj w razie potrzeby)

% Obliczenie rozmiarów okien
window_size = round(window_duration * fs);
step_size = round(step_duration * fs);
num_windows = floor((length(sx) - window_size) / step_size) + 1;

% Inicjalizacja wyników
detected_digits = [];

for i = 1:num_windows
    start_idx = (i-1) * step_size + 1;
    end_idx = start_idx + window_size - 1;
    if end_idx > length(sx)
        break;
    end
    frame = sx(start_idx:end_idx);  % <-- nowa nazwa zamiast "window"

    % Analiza przy użyciu filtrów IIR z jednym biegunem (Butterworth 1 rzędu)
    power_low = zeros(size(dtmf_freqs_low));
    power_high = zeros(size(dtmf_freqs_high));

    for k = 1:length(dtmf_freqs_low)
        fc = dtmf_freqs_low(k);
        [b, a] = butter(1, [fc-10 fc+10]/(fs/2), 'bandpass');  % (1 rzad - 1 biegun) wąskie pasmo ±10Hz
        filtered = filter(b, a, frame);
        power_low(k) = sum(filtered.^2);
    end

    for k = 1:length(dtmf_freqs_high)
        fc = dtmf_freqs_high(k);
        [b, a] = butter(1, [fc-10 fc+10]/(fs/2), 'bandpass');
        filtered = filter(b, a, frame);
        power_high(k) = sum(filtered.^2);
    end

    [max_power_low, idx_low] = max(power_low);
    [max_power_high, idx_high] = max(power_high);

    % Wybieramy częstotliwości o najwyższej energii
    low_freq = dtmf_freqs_low(idx_low);
    high_freq = dtmf_freqs_high(idx_high);

    % Przypisanie cyfry
    digit = map_dtmf(low_freq, high_freq);
    detected_digits = [detected_digits, digit];
end

% Wyświetlenie wykrytej sekwencji
disp(['Wykryta sekwencja z IIR dla ', fname, ': ', detected_digits]);



%% 9. Opcjonalne (+0,25) – Algorytm decyzyjny
% --- Parametry obwiedni ---
env_win_sec   = 0.02;   % [s] okno ruchomej średniej (20 ms)
min_tone_sec  = 0.15;    % minimalna długość tonu [s]
min_pause_sec = 0.2;    % minimalna przerwa między tonami [s]

env_win    = round(env_win_sec   * fs);
min_tone   = round(min_tone_sec  * fs);
min_pause  = round(min_pause_sec * fs);

% RMS-owa obwiednia
env = movmean(abs(sx), env_win);

% --- Maska i wycinanie segmentów ---
thr_env = 0.18 * max(env);        % próg (20% maks.) 
mask    = env > thr_env;
dmask   = diff([0; mask; 0]);

starts  = find(dmask==1);
ends    = find(dmask==-1)-1;

% usuń zbyt krótkie albo bez pauzy od poprzedniego
good = false(size(starts));
for k = 1:numel(starts)
    L = ends(k)-starts(k)+1;
    if L>=min_tone && (k==1 || starts(k)-ends(k-1)>=min_pause)
        good(k) = true;
    end
end
starts = starts(good);
ends   = ends(good);

% ogranicz do dokładnie 5 segmentów
n = min(5,numel(starts));
starts = starts(1:n);
ends   = ends(1:n);

% oblicz centra segmentów
centers = round((starts + ends)/2);

% --- (opcjonalnie) podgląd wykrycia ---
figure; plot((0:length(env)-1)/fs, env); hold on;
for k=1:n
    plot([starts(k), ends(k)]/fs, [thr_env thr_env],'r','LineWidth',2);
    plot(centers(k)/fs, env(centers(k)),'ko','MarkerFaceColor','w');
end
xlabel('t [s]'); ylabel('obwiednia'); title('Wykryte segmenty tonów');

% --- Goertzel na każdym środkowym punkcie ---
dtmf_low   = [697,770,852,941];
dtmf_high  = [1209,1336,1477];
thr_goer   = 0.5;           % próg mocy Goertzla
win_dur    = 0.6;           % długość okna [s]
win_size   = round(win_dur * fs);
half_size  = floor(win_size/2);

decoded = '';
for k = 1:n
    c  = centers(k);
    st = max(1, c-half_size);
    en = min(length(sx), st+win_size-1);
    w  = sx(st:en);

    pl = abs(goertzel(w, dtmf_low  * win_size/fs));
    ph = abs(goertzel(w, dtmf_high * win_size/fs));
    [ml, il] = max(pl);
    [mh, ih] = max(ph);

    if ml>thr_goer && mh>thr_goer
        decoded(end+1) = map_dtmf(dtmf_low(il), dtmf_high(ih));
    else
        decoded(end+1) = '?';
    end
end

disp(['Wykryta sekwencja: ', decoded]);

function d = map_dtmf(low, high)
    keys = {'697_1209','697_1336','697_1477','770_1209','770_1336', ...
            '770_1477','852_1209','852_1336','852_1477','941_1209', ...
            '941_1336','941_1477'};
    digs = {'1','2','3','4','5','6','7','8','9','*','0','#'};
    key  = sprintf('%d_%d', low, high);
    idx  = find(strcmp(keys, key),1);
    d    = digs{idx};
end

% -------------------------------------------------------------------------
% Czy analiza Goertzlem jest łatwiejsza?
% – nie wymaga filtrów BP ani pełnych spektrogramów (mniej pamięci i czas JIT),
% – pozwala łatwo „wycelować” dokładnie w interesujące nas częstotliwości,
% – wystarczy jedynie parę mnożeń/dodawań na próbkę → nadaje się do DSP real-time.
% -------------------------------------------------------------------------
