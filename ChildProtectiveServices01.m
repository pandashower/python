clear all; close all; clc;

%% Zad.1 (2,5p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% A)
t = 0.1; % s
A = 230; % V
f = 50; % Hz

fs1 = 10000;  % 10 kHz (pseudo-analog)
fs2 = 500;    % 500 Hz
fs3 = 200;    % 200 Hz

t1 = 0 : 1/fs1 : t; % 1/fs bo odlegość czasowa między kolejnymi próbkami
y1 = A * sin(2*pi*f*t1);

t2 = 0 : 1/fs2 : t;
y2 = A * sin(2*pi*f*t2);

t3 = 0 : 1/fs3 : t;
y3 = A * sin(2*pi*f*t3);

figure;
plot(t1, y1, 'b-');
hold on;
plot(t2, y2, 'r-o');
plot(t3, y3, 'k-x');
hold off;

xlabel('Czas [s]');
ylabel('Napięcie [V]');
title('Porównanie próbkowania sinusoidy o f=50 Hz przy różnych fs');
legend('f_s = 10 kHz','f_s = 500 Hz','f_s = 200 Hz','Location','best');
grid on;

% B)
t_end = 1; fs1 = 10000; fs2 = 51; fs3 = 50; fs4 = 49;
% fs2 = 26; fs3 = 25; fs4 = 24; % ODKOMENTOWAĆ DLA 2 WERSJI

t_prim1 = 0: 1/fs1 : t_end;
y_prim1 = A * sin(2*pi*f*t_prim1);

t_prim2 = 0 : 1/fs2 : t_end;
y_prim2 = A * sin(2*pi*f*t_prim2);

t_prim3 = 0 : 1/fs3 : t_end;
y_prim3 = A * sin(2*pi*f*t_prim3);

t_prim4 = 0 : 1/fs4 : t_end;
y_prim4 = A * sin(2*pi*f*t_prim4);

figure;
plot(t_prim1, y_prim1, 'b-');
hold on;
plot(t_prim2, y_prim2, 'g-o');
plot(t_prim3, y_prim3, 'r-o');
plot(t_prim4, y_prim4, 'k-o');
hold off;

xlabel('Czas [s]');
ylabel('Napięcie [V]');
title('Porównanie próbkowania sinusoidy o f=50 Hz przy różnych fs');
legend('f_s = 10 kHz','f_s = 51 Hz','f_s = 50 Hz','f_s = 49 Hz','Location','best');
grid on;


%% C)
clear all; close all; clc;

fs = 100;             % [Hz]
t_end = 1;            % [s]
t = 0 : 1/fs : t_end;

f_start = 0;          % [Hz]
f_end   = 300;         % [Hz]
f_step  = 5;        % [Hz]
freq_values = f_start : f_step : f_end;

for i = 1 : length(freq_values)
    f_sig = freq_values(i);                      % f sygnału w danej iteracji
    y = sin(2*pi*f_sig*t);                       % generujemy sinusoidę
    plot(t, y, 'b-');                            % rysujemy sinusoidę
    xlabel('Czas [s]');
    ylabel('Amplituda');

    disp(['Iteracja: ' num2str(i) ...
          ', Częstotliwość = ' num2str(f_sig) ' Hz']);

    pause(0.1);   % żeby było widać animację
end

% Porównania:
% czest Nyquista = 50Hz wszytkie powyżej aliasowane bo fs = 100Hz
f_test1 = [5, 105, 205]; % przez alising pkrywaja sie! f_alias = |f_sig-k*fs|, gdzie
% k to całkowita, która sprawia, że wynik jest od 0 do fs/2
f_test2 = [95, 195, 295]; % Zmiana znaku jest skutkiem tego, że częstotliwości
% aliasowane zachowują się jak (po powyzej f Nyqusita):
% sin(2π(f_s - f)t) = -sin(2π f t) 
f_test3 = [95, 105]; % podobnue jak wyżej

figure;
hold on;
plot(t, sin(2*pi*f_test1(1)*t), 'r',  'DisplayName','5 Hz');
plot(t, sin(2*pi*f_test1(2)*t), 'g',  'DisplayName','10 Hz');
plot(t, sin(2*pi*f_test1(3)*t), 'b--','DisplayName','105 Hz');
hold off;
xlabel('Czas [s]');
ylabel('Amplituda');
title('Porownanie sinusoid: 5 Hz, 105 Hz, 205 Hz (fs=100 Hz)');
legend('Location','best');
grid on;

figure;
hold on;
plot(t, sin(2*pi*f_test2(1)*t), 'r',  'DisplayName','95 Hz');
plot(t, sin(2*pi*f_test2(2)*t), 'g',  'DisplayName','195 Hz');
plot(t, sin(2*pi*f_test2(3)*t), 'b--','DisplayName','295 Hz');
hold off;
xlabel('Czas [s]');
ylabel('Amplituda');
title('Porownanie sinusoid: 95 Hz, 195 Hz, 295 Hz (fs=100 Hz)');
legend('Location','best');
grid on;

figure;
hold on;
plot(t, sin(2*pi*f_test3(1)*t), 'r',  'DisplayName','95 Hz');
plot(t, sin(2*pi*f_test3(2)*t), 'g',  'DisplayName','105 Hz');
hold off;
xlabel('Czas [s]');
ylabel('Amplituda');
title('Porownanie sinusoid: 95 Hz, 105 Hz(fs=100 Hz)');
legend('Location','best');
grid on;

% Dla cos:
figure;
hold on;
plot(t, cos(2*pi*f_test1(1)*t), 'r',  'DisplayName','5 Hz');
plot(t, cos(2*pi*f_test1(2)*t), 'g',  'DisplayName','10 Hz');
plot(t, cos(2*pi*f_test1(3)*t), 'b--','DisplayName','105 Hz');
hold off;
xlabel('Czas [s]');
ylabel('Amplituda');
title('Porownanie cosinusoid: 5 Hz, 105 Hz, 205 Hz (fs=100 Hz)');
legend('Location','best');
grid on;

figure;
hold on;
plot(t, cos(2*pi*f_test2(1)*t), 'r',  'DisplayName','95 Hz');
plot(t, cos(2*pi*f_test2(2)*t), 'g',  'DisplayName','195 Hz');
plot(t, cos(2*pi*f_test2(3)*t), 'b--','DisplayName','295 Hz');
hold off;
xlabel('Czas [s]');
ylabel('Amplituda');
title('Porownanie cosinusoid: 95 Hz, 195 Hz, 295 Hz (fs=100 Hz)');
legend('Location','best');
grid on;

figure;
hold on;
plot(t, cos(2*pi*f_test3(1)*t), 'r',  'DisplayName','95 Hz');
plot(t, cos(2*pi*f_test3(2)*t), 'g',  'DisplayName','105 Hz');
hold off;
xlabel('Czas [s]');
ylabel('Amplituda');
title('Porownanie cosinusoid: 95 Hz, 105 Hz(fs=100 Hz)');
legend('Location','best');
grid on;


%% D)
clear all; close all; clc;
fs = 10000; % częstotliwość próbkowania (10 kHz)
t_end = 1; % czas trwania sygnału [s]
t = 0 : 1/fs : t_end; % oś czasu 
fn = 50; % częstotliwość nośna [Hz]
fm = 1; % częstotliwość modulująca [Hz]
dev = 5; % głębokość modulacji - tu dewiacja raczej (±5 Hz)
A = 1; % amplituda sygnału - zakładam 1 bo brak w zad.

% Sygnał modulujący (poglądowo)
m_t = sin(2*pi*fm*t);

% Sygnał z modulacją częstotliwości (SFM) - zmodulowany:
beta  = dev / fm; % współczynnik modulacji (głębokość)
y_mod = A * sin( 2*pi*fn*t + beta*sin(2*pi*fm*t) );

figure('Name','Sygnał zmodulowany i modulujący');
subplot(2,1,1);
plot(t, m_t, 'r'); 
xlabel('Czas [s]'); ylabel('Amplituda'); 
title('Sygnał modulujący (1 Hz)');
grid on;
subplot(2,1,2);
plot(t, y_mod, 'b'); 
xlabel('Czas [s]'); ylabel('Amplituda');
title('Sygnał zmodulowany (SFM, nośna 50 Hz, dev ±5 Hz)');
grid on;

% [2] Próbkowanie z fs2 = 25 Hz i porównanie
fs2 = 25;
t2  = 0 : 1/fs2 : t_end;
y_samp = A * sin( 2*pi*fn*t2 + beta*sin(2*pi*fm*t2) );

figure('Name','Porównanie: sygnał oryginalny i próbkowany');
subplot(3,1,1);
plot(t, y_mod, 'DisplayName','Oryginał (fs=10 kHz)');
hold on;
stem(t2, y_samp, 'r', 'DisplayName','Próbki (fs2=25 Hz)');
hold off;
legend('Location','best'); grid on;
xlabel('Czas [s]'); ylabel('Amplituda');
title('Sygnał zmodulowany i jego próbki');

% Rekonstrukcja sygnału (Zero-Order Hold - ZOH)
y_rec = interp1(t2, y_samp, t, 'previous', 'extrap');
err = y_mod - y_rec;

subplot(3,1,2);
plot(t, y_mod, 'b', 'DisplayName','Oryginał');
hold on;
plot(t, y_rec, 'm--', 'DisplayName','Rekonstrukcja (ZOH)');
hold off;
legend('Location','best'); grid on;
xlabel('Czas [s]'); ylabel('Amplituda');
title('Porównanie oryginału z rekonstrukcją');

subplot(3,1,3);
plot(t, err);
grid on;
xlabel('Czas [s]'); ylabel('Błąd');
title('Błąd odtwarzania (oryginał - rekonstrukcja)');

% [3] Widmo gęstości mocy (PSD) przed i po próbkowaniu
figure;
subplot(2,1,1);
pwelch(y_mod, [], [], [], fs);
title('Widmo gęstości mocy - sygnał przed próbkowaniem');

subplot(2,1,2);
pwelch(y_samp, [], [], [], fs2);
title('Widmo gęstości mocy - sygnał po próbkowaniu');




%% Zad.2 (2p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% Zad.1 A)
t = 0.1; % s
A = 230; % V
f = 50; % Hz
fs = 10000;  % 10 kHz (pseudo-analog)
fs2 = 500;    % 500 Hz
fs3 = 200;    % 200 Hz
t1 = 0 : 1/fs : t; % 1/fs bo odlegość czasowa między kolejnymi próbkami
y1 = A * sin(2*pi*f*t1);
t2 = 0 : 1/fs2 : t;
y2 = A * sin(2*pi*f*t2);
t3 = 0 : 1/fs3 : t;
y3 = A * sin(2*pi*f*t3);
% KONIEC ZAD.1

% Rekonstrukcja sygnału za pomocą splotu z sinc:
% wychodzimy od sygnału y3 próbkowanego z fs3 = 200 Hz (t3, y3)
% chcemy go odtworzyć na gęstszej siatce czasowej (t1) używanej w "pseudo-analogu" fs=10 kHz.

% okres próbkowania dla sygnału o fs3 = 200 Hz
T3 = 1/fs3;

% wybieramy siatkę czasową, t1 żeby móc łatwo porównać z "pseudo-analogiem".
ts = t1;  

% wektor na zrekonstruowany sygnał:
x_hat = zeros(size(ts));

% Główna pętla po chwilach, w których chcemy obliczyć wartość zrekonstruowaną
for i = 1:length(ts)
    % aktualny czas, w którym odtwarzamy sygnał
    ti = ts(i);
    
    % sumujemy wkłady od wszystkich próbek y3(n) * sinc( (t - nT)/T )
    for n = 1:length(t3)
        % (ti - t3(n)) / T3 -> argument dla funkcji sinc
        x_hat(i) = x_hat(i) + y3(n) * sinc( (ti - t3(n)) / T3 );
    end
end

figure('Name','Rekonstrukcja sygnału z fs3=200Hz za pomocą sinc');
plot(t1, y1, 'b', 'LineWidth',1.5);  hold on;
plot(ts, x_hat, 'r--', 'LineWidth',1.2);
xlabel('Czas [s]'); ylabel('Amplituda');
legend('pseudo-analog fs=10kHz','rekonstrukcja sinc (fs3=200Hz)');
title('Rekonstrukcja metodą splotu z sinc');

blad = x_hat - y1;  % różnica między zrekonstruowanym a "pseudo-analogiem"
figure('Name','Błąd rekonstrukcji');
plot(ts, blad);
xlabel('Czas [s]'); ylabel('Błąd');
title('Błąd rekonstrukcji (sygnał odtworzony - pseudo-analog)');
grid on;




%% Zad.3 (1,5p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
function [r, lags] = my_xcorr(x, y)
    Nx = length(x);
    Ny = length(y);

    % Długość korelacji
    r = zeros(1, Nx + Ny - 1);
    lags = -(Ny-1):(Nx-1); % Obliczanie lagów

    for lag = -(Ny-1):(Nx-1)
        suma = 0;
        for i = 1:Ny
            j = i + lag;
            if j > 0 && j <= Nx
                suma = suma + x(j) * y(i);
            end
        end
        r(lag + Ny) = suma;
    end
end

load('CPS01\adsl_x.mat'); % w pliku jest zmienna 'x'

M = 32;     % długość prefiksu
N = 512;    % długość ramki

% 2. Definiujemy wzorzec: ostatnie M próbek z ramki
%   (przyjmujemy, że ramka startuje od x(1) do x(N))
prefix_ref = x(N - M + 1 : N);

% 3. Obliczamy korelację krzyżową sygnału x z wzorcem prefix_ref
[R, lags] = xcorr(x, prefix_ref);

% 4. Wyszukujemy lokalne maksima korelacji (pozwala to znaleźć wystąpienia prefiksu)
[pkVals, pkLocs] = findpeaks(R, ...
    'MinPeakHeight', 0.5*max(R)); 
%  - 'MinPeakHeight' można dostosować zależnie od poziomu szumu/energii w sygnale

% Indeksy w 'lags' odpowiadają przesunięciom między x a prefix_ref.
% Żeby uzyskać indeks w sygnale x, dodajemy długość prefiksu (lub N, w zależności
% od interpretacji). Najczęściej wystąpienie prefiksu zaczyna się przy:
startIndices = lags(pkLocs) + length(prefix_ref);

% Wyświetlamy wyniki
disp('Znalezione indeksy początków prefiksów w sygnale:');
disp(startIndices);

% (Opcjonalnie) Rysujemy korelację i zaznaczamy piki
figure; 
plot(lags, R); hold on;
plot(lags(pkLocs), pkVals, 'ro', 'MarkerFaceColor','r');
xlabel('Lag'); ylabel('Wartość korelacji');
title('Korelacja krzyżowa sygnału z wzorcem prefiksu');
grid on;




%% Zad.4 (1p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
f_pr = 16000;  % Częstotliwość próbkowania w Hz (16 kHz)
T = 0.1;       % Czas trwania bitu w sekundach (100 ms)
f_c = 500;     % Częstotliwość nośna w Hz
name = 'Mateusz';

% Konwersja imienia na kody ASCII
ascii_codes = double(name);

% Konwersja kodów ASCII na binarne (8 bitów na znak)
binary_str = dec2bin(ascii_codes, 8);
bits = binary_str(:)' - '0';  % Konwersja na tablicę numeryczną bitów

% Obliczenie liczby próbek na bit
N = round(f_pr * T);  % Liczba próbek na bit

% Wektor czasu dla jednego bitu
t = (0:N-1) / f_pr;

% Inicjalizacja sygnału
signal = [];

% Generowanie sygnału dla każdego bitu
for bit = bits
    if bit == 0
        segment = sin(2 * pi * f_c * t);
    else
        segment = -sin(2 * pi * f_c * t);
    end
    signal = [signal, segment];
end

t_total = (0:length(signal)-1) / f_pr;

% Narysowanie wykresu sygnału
figure; % Nowe okno dla wykresu
plot(t_total, signal);
xlabel('Czas [s]');
ylabel('Amplituda');
title('Sygnał transmitujący bity kodu ASCII imienia');
grid on; % Włączenie siatki dla lepszej czytelności

% Definiowanie częstotliwości próbkowania do odtwarzania
sampling_freqs = [8000, 16000, 24000, 32000, 48000];

% Odtwarzanie sygnału przy różnych częstotliwościach
for fs = sampling_freqs
    disp(['Odtwarzanie przy częstotliwości próbkowania: ', num2str(fs), ' Hz']);
    soundsc(signal, fs);
    pause(length(signal)/fs + 1);  % Pauza, aby sygnał mógł się w pełni odtworzyć
end
