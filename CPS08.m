%%  Lab 08 — Hilbert FIR + demodulacja AM
%  Wersja „windowed‑sinc + okno Blackmana” w pełni zgodna z opisem rys. 12.18
%  -------------------------------------------------------------------------
%  • Filtr Hilberta budowany z idealnej odpowiedzi 2/(π n) przyciętej i
%    pomnożonej przez okno Blackmana.
%  • Zachowano sym‑padding, korektę |H|, LS‑dopasowanie amplitud.
%  • Wymuszamy wektory kolumnowe, aby uniknąć błędu z vertcat.
%  -------------------------------------------------------------------------

clear;  clc;  close all

%% 0)  Parametry podstawowe
fs = 1e3;       % [Hz] — częstotliwość próbkowania
fc = 200;       % [Hz] — częstotliwość nośna
T  = 1;         % [s ] — czas trwania sygnału
N  = fs*T;      % liczba próbek

t  = (0:N-1).'/fs;              % wektor czasu (kolumna)

%% 1)  Wczytaj sygnał i upewnij się, że to kolumna
load lab08_am.mat               % plik z wieloma realizacjami (s1…sX)
idx = 4;                         % <— wybierz numer realizacji
x   = eval(sprintf('s%d',idx));
x   = x(:);                      % wektor kolumnowy → spójny wymiar

%% 2)  Projekt FIR‑Hilbert metodą windowed‑sinc + Blackman
M   = 127;                       % długość filtru (liczba współczynników)
mid = (M-1)/2;                   % indeks n = 0 po przesunięciu
n   = (-mid:mid).';              % kolumna z próbkami n

h_ideal      = zeros(M,1);
oddSamples   = mod(n,2) ~= 0;     % jedynie nieparzyste n
h_ideal(oddSamples) = 2./(pi*n(oddSamples));  % 2/(π n) dla n nieparzystych
% n parzyste i 0 pozostają 0 —> przesunięcie fazy -90°

win = blackman(M);               % okno Blackmana (kolumna)
h   = h_ideal .* win;            % odpowiedź impulsowa filtru

grpDelay = mid;                  % opóźnienie grupowe = (M-1)/2

%  (istotne przy kompensacji fazy / obwiedni)
[H,w] = freqz(h,1,4096,fs);
Hc    = abs(interp1(w,H,fc));     % |H(jΩ_c)| ~ 1 (kontrola wzmocnienia)

%% 3)  Transformacja Hilberta z sym‑paddingiem
xPad = [ flipud(x(1:grpDelay));    % lewy padding (odbicie)
         x;                        % oryginalne próbki
         flipud(x(end-grpDelay+1:end)) ];  % prawy padding

xHilPad = filter(h,1,xPad);       % filtracja FIR
xHilPad = circshift(xHilPad,-grpDelay);    % kompensacja opóźnienia
xHil    = xHilPad(grpDelay+1 : grpDelay+N);% odcięcie padd.

amp = sqrt(x.^2 + xHil.^2) / Hc;  % obwiednia + kompensacja wzm.

%% 4)  Oddziel składową stałą
amp_dc = mean(amp);
m       = amp - amp_dc;           % sygnał modulujący (bez DC)

%% 5)  Widmo obwiedni → częstotliwości + amplitudy
Lfft   = 2^nextpow2(N*4);         % zero‑padding (x4) → Δf ≈ 0.24 Hz
winFFT = blackman(N);
Mwin   = sum(winFFT);

SPEC   = fft(m .* winFFT, Lfft);
f      = (0:Lfft-1).' * fs/Lfft;

mask = f < (fc-5);                % pasmo < fc − 5 Hz
[~,locs] = findpeaks(abs(SPEC(mask)), 'NPeaks',3,'SortStr','descend');
locs      = locs(:);
f_est     = f(locs);
[f_est,ix]= sort(f_est);          % rosnąco
locs      = locs(ix);

A_est = 2*abs(SPEC(locs)) / Mwin; % korekta okna Blackmana

% LS‑refinement (dokładniejsze amplitudy)
C = [cos(2*pi*f_est(1)*t), cos(2*pi*f_est(2)*t), cos(2*pi*f_est(3)*t)];
theta = C \ m;                     % least‑squares
A_est = theta(:).';

%% 6)  Rekonstrukcja i MSE
m_rec = amp_dc + A_est(1)*cos(2*pi*f_est(1)*t) + ...
                      A_est(2)*cos(2*pi*f_est(2)*t) + ...
                      A_est(3)*cos(2*pi*f_est(3)*t);

x_rec = m_rec .* cos(2*pi*fc*t);
MSE   = mean((x - x_rec).^2);

%% 7)  Raport
fprintf('\n=== Parametry modulacji (Hilbert + Blackman windowed‑sinc) ===\n');
for k = 1:3
    fprintf('   f%d = %6.2f Hz   A%d = %.4f\n', k, f_est(k), k, A_est(k));
end
fprintf('   skł. stała = %.4f\n', amp_dc);
fprintf('   MSE rekonstrukcji: %.3e\n', MSE);

%% 8)  Wykresy (opcjonalne)
figure;
subplot(2,1,1), plot(t, amp), grid on, title('Obwiednia')
subplot(2,1,2), plot(t,x,'b', t,x_rec,'r--'), grid on
legend('x – oryginał','x_{rec} – rekonstrukcja'), xlabel('t [s]')


















clear all; close all; clc;

%% Parametry
fs = 400e3; % częstotliwość próbkowania sygnału radiowego
fc1 = 100e3; % nośna 1
fc2 = 110e3; % nośna 2
dA = 0.25; % głębokość modulacji

%% Wczytanie pliku audio, nroma i nadpróbkowanie
[x1, fsx] = audioread('mowa8000.wav');
x2 = flipud(x1); % odwrotnie puszczona mowa

% Normalizacja
x1 = x1 / max(abs(x1));
x2 = x2 / max(abs(x2));

% Nadpróbkowanie
L = fs / fsx;
x1u = resample(x1, fs, fsx);
x2u = resample(x2, fs, fsx);

t = (0:length(x1u)-1)'/fs;

%% MODULACJA DSB-C
y1 = (1 + dA * x1u) .* cos(2*pi*fc1*t);
y2 = (1 + dA * x2u) .* cos(2*pi*fc2*t);
yDSBC = y1 + y2;

%% MODULACJA DSB-SC
y1 = dA * x1u .* cos(2*pi*fc1*t);
y2 = dA * x2u .* cos(2*pi*fc2*t);
yDSBSC = y1 + y2;

%% MODULACJA SSB-SC
% Filtr Hilberta FIR (np. 101 współczynników)
N = 101;
h = firpm(N-1, [0.05 0.95], [1 1], 'hilbert'); % filtr pasmowo-przepustowy

x1H = filter(h, 1, x1u);
x2H = filter(h, 1, x2u);

% SSB-SC: jedna prawa (USB), jedna lewa (LSB)
y1 = 0.5 * dA * x1u .* cos(2*pi*fc1*t) - 0.5 * dA * x1H .* sin(2*pi*fc1*t); % LSB
y2 = 0.5 * dA * x2u .* cos(2*pi*fc2*t) + 0.5 * dA * x2H .* sin(2*pi*fc2*t); % USB
ySSBSC = y1 + y2;


%% Projektowanie filtru Hilberta (FIR)
N = 201;                   % długość filtru (nieparzysta liczba próbek, typowe dla filtru FIR)
n = -(N-1)/2:(N-1)/2;     % oś próbkowania odpowiedzi impulsowej (indeksy próbek)

% Tworzenie idealnej odpowiedzi impulsowej filtru Hilberta
h_ideal = (1 ./ (pi * n)) .* (1 - cos(pi * n));  % h[n] = 1/(pi*n) * (1 - cos(pi*n))
h_ideal((N+1)/2) = 0;     % ustawienie wartości w n=0 (uniknięcie NaN) zgodnie z definicją

% Zastosowanie okna Blackmana (poprawa charakterystyki filtru)
h = h_ideal .* blackman(N)';  % Okno 

% Normalizacja filtru Hilberta
h = h / sum(h);  % Normalizacja, aby suma współczynników była równa 1 (zapewnia to stabilność)

%% PRZYGOTOWANIE – DEMODULACJA dla każdej transmisji z użyciem filtru Hilberta FIR
hilbert_filt = @(sig) real(conv(sig, h, 'same')); % Filtr Hilberta przy użyciu splotu

% DSB-C
demod_DSB_C1 = resample(hilbert_filt(yDSBC .* cos(2*pi*fc1*t)), fsx, fs);
demod_DSB_C2 = resample(hilbert_filt(yDSBC .* cos(2*pi*fc2*t)), fsx, fs);
demod_DSB_C2 = flipud(demod_DSB_C2); % cofamy odwrócenie

% DSB-SC
demod_DSB_SC1 = resample(hilbert_filt(yDSBSC .* cos(2*pi*fc1*t)), fsx, fs);
demod_DSB_SC2 = resample(hilbert_filt(yDSBSC .* cos(2*pi*fc2*t)), fsx, fs);
demod_DSB_SC2 = flipud(demod_DSB_SC2);

% SSB-SC
ssb_demod1 = ySSBSC .* (cos(2*pi*fc1*t) - 1i*sin(2*pi*fc1*t)); % LSB (x1)
ssb_demod2 = ySSBSC .* (cos(2*pi*fc2*t) + 1i*sin(2*pi*fc2*t)); % USB (x2)
demod_SSB_SC1 = resample(hilbert_filt(ssb_demod1), fsx, fs);
demod_SSB_SC2 = resample(hilbert_filt(ssb_demod2), fsx, fs);
demod_SSB_SC2 = flipud(demod_SSB_SC2);

%% ODTWARZANIE
fprintf("\nOdtwarzanie transmisji DSB-C...\n");
disp("Stacja 1 – po demodulacji");
soundsc(demod_DSB_C1, fsx); pause();

disp("Stacja 2 – po demodulacji");
soundsc(demod_DSB_C2, fsx); pause();

fprintf("\nOdtwarzanie transmisji DSB-SC...\n");
disp("Stacja 1 – po demodulacji");
soundsc(demod_DSB_SC1, fsx); pause();

disp("Stacja 2 – po demodulacji");
soundsc(demod_DSB_SC2, fsx); pause();

fprintf("\nOdtwarzanie transmisji SSB-SC...\n");
disp("Stacja 1 – po demodulacji");
soundsc(demod_SSB_SC1, fsx); pause();

disp("Stacja 2 – po demodulacji");
soundsc(demod_SSB_SC2, fsx); pause();


%% ----------------------------------------------------------------------
%% 4.  OPCJA 2 – dwie stacje na wspólnej nośnej (SSB)  – poprawka
%% ----------------------------------------------------------------------
lowpass4k = @(s) lowpass(s, 4e3, fs);       % dolnoprzepustowy 4 kHz
analytic = @(sig) sig + 1j*hilbert_filt(sig);   % używa mojego filtru h
% ---- MODULACJA ----
yLSB = 0.5*x1u .* cos(2*pi*fc1*t) - 0.5*x1H .* sin(2*pi*fc1*t);
yUSB = 0.5*x2u .* cos(2*pi*fc1*t) + 0.5*x2H .* sin(2*pi*fc1*t);
ySSBSC_dual = yLSB + yUSB;

% demodulacja
z = analytic(ySSBSC_dual) .* exp(-1j*2*pi*fc1*t);   % zejście do basebandu
bb_usb = lowpass4k(  real( z) );                    % +4 kHz  → x2  (USB)
bb_lsb = lowpass4k( -imag( z) );                    % –4 kHz  → x1  (LSB)

dec_usb = bb_usb  / max(abs(bb_usb));
dec_lsb = bb_lsb  / max(abs(bb_lsb));

dec_usb_down = resample( dec_usb, fsx, fs);
dec_lsb_down = resample( dec_lsb, fsx, fs);  % odwracamy przy odtwarzaniu

% Odsłuch
gain = 0.9;
fprintf('USB...\n');  sound(dec_usb_down * gain, fsx);
pause;
fprintf('LSB (odwrócona mowa)...\n');
sound(flipud(dec_lsb_down) * gain, fsx);













clear all; close all; clc;

%% 1. Parametry sygnałów
fs_target = 48000; % Docelowa częstotliwość próbkowania
t = 1; % Czas trwania sygnału

% Częstotliwości i próbkowania sygnałów
f1 = 1001.2; fs1 = 8000;
f2 = 303.1;  fs2 = 32000;
f3 = 2110.4; fs3 = 48000;

%% 2. Generacja sygnałów sinusoidalnych
t1 = (0:1/fs1:t-1/fs1)'; x1 = sin(2*pi*f1*t1);
t2 = (0:1/fs2:t-1/fs2)'; x2 = sin(2*pi*f2*t2);
t3 = (0:1/fs3:t-1/fs3)'; x3 = sin(2*pi*f3*t3);
t_target = (0:1/fs_target:t-1/fs_target)';
x4_expected = sin(2*pi*f1*t_target) + sin(2*pi*f2*t_target) + sin(2*pi*f3*t_target);

%% 3. Repróbkowanie metodą resample
[p1, q1] = rat(fs_target/fs1);
[p2, q2] = rat(fs_target/fs2);
x1_resampled = resample(x1, p1, q1);
x2_resampled = resample(x2, p2, q2);
x3_resampled = x3;

% Dopasowanie długości
min_len = min([length(x1_resampled), length(x2_resampled), length(x3_resampled)]);
x1_resampled = x1_resampled(1:min_len);
x2_resampled = x2_resampled(1:min_len);
x3_resampled = x3_resampled(1:min_len);
x4_resampled = x1_resampled + x2_resampled + x3_resampled;
x4_expected = x4_expected(1:min_len);

%% 4. Obliczenie błędu MSE
mse = mean((x4_resampled - x4_expected).^2);
fprintf('MSE między sygnałem wynikowym a oczekiwanym: %.10e\n', mse);

%% 5. Wykresy porównawcze
samples_to_show = 500;

figure('Name', 'Porównanie repróbkowanych sygnałów', 'Position', [100, 100, 900, 700]);
subplot(3,1,1);
plot(t_target(1:samples_to_show), x1_resampled(1:samples_to_show), 'r');
title('Sygnał x1 po resample'); grid on;

subplot(3,1,2);
plot(t_target(1:samples_to_show), x2_resampled(1:samples_to_show), 'g');
title('Sygnał x2 po resample'); grid on;

subplot(3,1,3);
plot(t_target(1:samples_to_show), x4_resampled(1:samples_to_show), 'b', ...
     t_target(1:samples_to_show), x4_expected(1:samples_to_show), 'k--');
title('Sygnał wynikowy vs. oczekiwany'); legend('Wynikowy', 'Oczekiwany'); grid on;

%% 6. Analiza widmowa
figure('Name', 'Widmo sygnału wynikowego', 'Position', [100, 100, 900, 500]);
NFFT = 2^nextpow2(length(x4_resampled));
X4 = fft(x4_resampled, NFFT) / length(x4_resampled);
f = fs_target/2 * linspace(0, 1, NFFT/2+1);
plot(f, 2*abs(X4(1:NFFT/2+1)));
title('Widmo sygnału x4'); xlabel('Hz'); ylabel('Amplituda'); grid on; xlim([0, 3000]);
hold on;
line([f1 f1], [0 0.5], 'Color', 'r', 'LineStyle', '--');
line([f2 f2], [0 0.5], 'Color', 'g', 'LineStyle', '--');
line([f3 f3], [0 0.5], 'Color', 'b', 'LineStyle', '--');

%% 7. Odsłuch
x4_expected = x4_expected / max(abs(x4_expected));
x4_play = x4_resampled / max(abs(x4_resampled));
fprintf('oczekiwany, ')
sound(x4_expected, fs_target); pause;
fprintf('otrzymany')
sound(x4_play, fs_target); pause;

%% 8. Miksowanie plików WAV  → 48 kHz 
try
    [wav1, fs_wav1] = audioread('x1.wav');
    [wav2, fs_wav2] = audioread('x2.wav');

    % --- MONO sumowanie kanałów stereo -------------------------------
    if size(wav1,2) > 1, wav1 = mean(wav1,2); end
    if size(wav2,2) > 1, wav2 = mean(wav2,2); end

    % ======================== 48 kHz =================================
    fs_target = 48e3;
    [p1,q1] = rat(fs_target/fs_wav1);
    [p2,q2] = rat(fs_target/fs_wav2);
    wav1_48 = resample(wav1,p1,q1);
    wav2_48 = resample(wav2,p2,q2);

    minLen = min(length(wav1_48),length(wav2_48));
    mix_48k = wav1_48(1:minLen) + wav2_48(1:minLen);
    mix_48k = mix_48k ./ max(abs(mix_48k));      % normalizacja
   
    % --- Odsłuch ----------------------------
    fprintf(' Odtwarzanie 48k\n');
    sound(mix_48k,  48e3);  

catch ME
    warning('Miksowanie przerwane: %s');
end
















%% FM – filtr różniczkujący z pasmem  (Zadanie 4)
%  autor: <Twoje imię>, data: <dzisiejsza>

clear all; clc; close all;

%% --- PARAMETRY ----------------------------------------------------------
fs       = 2e6;      % próbkowanie (2 MHz)
fc       = 200e3;    % „nowa” nośna po podpróbkowaniu (200 kHz)
% pasmo stacji FM (przyjęte roboczo):
pb1      = 150e3;    % dolna krawędź pasma użytecznego
pb2      = 250e3;    % górna krawędź pasma użytecznego
sb1      = 120e3;    % dolna krawędź zaporowa
sb2      = 280e3;    % górna krawędź zaporowa
astop    = 80;       % tłumienie w zaporze [dB]

audioFs  = 8e3;      % docelowe fs dla sygnału mowy
%% ------------------------------------------------------------------------

%% 1. Wczytanie danych -----------------------------------------------------
load lab08_fm.mat    % w pliku jest wektor x_fm (kompleksowy lub rzeczywisty)
x = x(:);         % upewniamy się, że to kolumna
N  = length(x);

%% 2A. Metoda z kaskadą dwóch filtrów  ------------------------------------
%  -- projekt DIFF (pełnopasmowy) --
Ndiff = 129;   % rząd (nieparzysty ➔ faza liniowa)
hdiff = firls(Ndiff-1,[0 1],[0 1],'differentiator');

%  -- projekt pasma BP --
edges = [0 sb1 pb1 pb2 sb2 fs/2] / (fs/2);   % normalizacja do Nyquista
amps  = [0  0   1   1   0   0];
Nbp   = 513;                                 % rząd zapewniający ~80 dB
hbp   = firls(Nbp-1,edges,amps);

%  -- filtr składający (splot) --
hcasc = conv(hdiff, hbp);

%% 2B. (1 pkt) Filtr różniczkujący pasmowo‑przepustowy w 1 kroku -----------
%   ①  krawędzie NIE zawierają 0 Hz
%   ②  amplitudy są proporcjonalne do F
%   ③  dodatkowe wagi mocno dociskają zaporę
%  Częstotliwości znormalizowane do Nyquista (fs/2)

pb1n = pb1/(fs/2);           % 0.15
pb2n = pb2/(fs/2);           % 0.25
sb1n = sb1/(fs/2);           % 0.12
sb2n = sb2/(fs/2);           % 0.28

%  F  – musi ZACZYNAĆ SIĘ od 0 i KOŃCZYĆ na 1
edgesDir = [0 sb1n  pb1n  pb2n  sb2n  1];

%  A  – dla trybu 'differentiator' podajemy |H(f)| / f
%       Chcemy |H(f)| = k·f w paśmie, więc  A = k = const = 1
ampsDir  = [0  0     1     1     0     0];

%  W  – wygórowane wagi w zaporach gwarantują ≥ 80 dB
weights  = [40  1  40];        % [stop  pass  stop]

Ndir = 1201;                   % rząd ~~ 80 dB (możesz zwiększyć do 1501)
hdir = firls(Ndir-1, edgesDir, ampsDir, weights, 'differentiator');

%% 3. Filtracja (obie metody do porównania) --------------------------------
y_casc = filter(hcasc,1,x);        % kaskada
y_dir  = filter(hdir ,1,x);        % pojedynczy filtr (1 pkt)

%% 4. Detekcja obwiedni, LPF, dekymacja ------------------------------------
%    (wykorzystujemy hermitowską reprezentację analityczną)
env_casc = abs(hilbert(y_casc));
env_dir  = abs(hilbert(y_dir ));

%  -- LPF do 15 kHz (≥ 16 kHz stop band) --
lpf = designfilt('lowpassfir', ...
                 'PassbandFrequency', 15e3, ...
                 'StopbandFrequency', 20e3, ...
                 'PassbandRipple', 1, ...
                 'StopbandAttenuation', 80, ...
                 'SampleRate', fs);
audio_casc = filtfilt(lpf, env_casc);   % zero‑fazowe, brak opóźnienia
audio_dir  = filtfilt(lpf, env_dir );

%  -- dekymacja do 8 kHz --
audio_casc = resample(audio_casc,audioFs,fs);
audio_dir  = resample(audio_dir ,audioFs,fs);

%% 5. Normalizacja głośności i zapis WAV ----------------------------------
audio_casc = audio_casc/max(abs(audio_casc))*0.99;
audio_dir  = audio_dir /max(abs(audio_dir ))*0.99;

audiowrite('mowa_fm_kaskada.wav',audio_casc,audioFs);
audiowrite('mowa_fm_firls.wav'   ,audio_dir ,audioFs);

fprintf('Gotowe! Pliki mowa_fm_*.wav zapisane w bieżącym katalogu.\n');

%% 6. Opcjonalna weryfikacja częstotliwościowa -----------------------------
%  (porównanie charakterystyk obu filtrów)
Nfft = 8192;
[Hcasc,w] = freqz(hcasc,1,Nfft,fs);
[Hdir ,~] = freqz(hdir ,1,Nfft,fs);

figure;
plot(w,20*log10(abs(Hcasc))+eps,'b'); hold on;
plot(w,20*log10(abs(Hdir ))+eps,'r--','LineWidth',1);
grid on; xlabel('f [Hz]'); ylabel('|H(f)| [dB]');
title('Charakterystyka różniczkującego BP: kaskada vs firls');
legend('DIFF + BP','firls(`differentiator`)','Location','Best');
xlim([0 500e3]); ylim([-140 20]);
