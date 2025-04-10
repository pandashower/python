clear all; close all; clc;

% --- Definicja specyfikacji filtra ---
fs = 256000; % Częstotliwość próbkowania [Hz]
fp = 64000;  % Górna częstotliwość pasma przepustowego (f_3dB) [Hz]
fsb = 128000; % Dolna częstotliwość pasma zaporowego (fs/2) [Hz]
Rp = 3;      % Maksymalne tłumienie (zafalowania) w paśmie przepustowym [dB]
Rs = 40;     % Minimalne tłumienie w paśmie zaporowym [dB]

% Konwersja częstotliwości do pulsacji [rad/s]
Wp = 2 * pi * fp; % Pulsacja graniczna pasma przepustowego
Ws = 2 * pi * fsb;% Pulsacja graniczna pasma zaporowego

% 1. Filtr Butterwortha
[n_butt, Wn_butt] = buttord(Wp, Ws, Rp, Rs, 's'); % Wyznacz rząd i częstotliwość graniczną
[z_butt, p_butt, k_butt] = butter(n_butt, Wn_butt, 's'); % Zaprojektuj filtr
sys_butt = zpk(z_butt, p_butt, k_butt); % Stwórz obiekt systemu ZPK
fprintf('Filtr Butterwortha: rząd n = %d\n', n_butt);

% 2. Filtr Czebyszewa typu 1
[n_cheb1, Wp_cheb1] = cheb1ord(Wp, Ws, Rp, Rs, 's'); % Wyznacz rząd i częstotliwość graniczną
[z_cheb1, p_cheb1, k_cheb1] = cheby1(n_cheb1, Rp, Wp_cheb1, 's'); % Zaprojektuj filtr
sys_cheb1 = zpk(z_cheb1, p_cheb1, k_cheb1); % Stwórz obiekt systemu ZPK
fprintf('Filtr Czebyszewa typu 1: rząd n = %d\n', n_cheb1);

% 3. Filtr Czebyszewa typu 2
[n_cheb2, Ws_cheb2] = cheb2ord(Wp, Ws, Rp, Rs, 's'); % Wyznacz rząd i częstotliwość graniczną
[z_cheb2, p_cheb2, k_cheb2] = cheby2(n_cheb2, Rs, Ws_cheb2, 's'); % Zaprojektuj filtr
sys_cheb2 = zpk(z_cheb2, p_cheb2, k_cheb2); % Stwórz obiekt systemu ZPK
fprintf('Filtr Czebyszewa typu 2: rząd n = %d\n', n_cheb2);

% 4. Filtr Eliptyczny (Cauer'a)
[n_ellip, Wp_ellip] = ellipord(Wp, Ws, Rp, Rs, 's'); % Wyznacz rząd i częstotliwość graniczną
[z_ellip, p_ellip, k_ellip] = ellip(n_ellip, Rp, Rs, Wp_ellip, 's'); % Zaprojektuj filtr
sys_ellip = zpk(z_ellip, p_ellip, k_ellip); % Stwórz obiekt systemu ZPK
fprintf('Filtr Eliptyczny: rząd n = %d\n', n_ellip);

% --- Rysowanie rozkładu biegunów i zer ---
figure;
subplot(2,2,1);
pzmap(sys_butt);
title(sprintf('Butterworth (n=%d)', n_butt));
grid on;

subplot(2,2,2);
pzmap(sys_cheb1);
title(sprintf('Czebyszew I (n=%d)', n_cheb1));
grid on;

subplot(2,2,3);
pzmap(sys_cheb2);
title(sprintf('Czebyszew II (n=%d)', n_cheb2));
grid on;

subplot(2,2,4);
pzmap(sys_ellip);
title(sprintf('Eliptyczny (n=%d)', n_ellip));
grid on;
sgtitle('Rozkład biegunów i zer zaprojektowanych filtrów'); % Tytuł główny

% --- Rysowanie charakterystyk częstotliwościowych ---
% Zakres częstotliwości do rysowania [Hz]
f_vec = linspace(0, fs, 1000); % Liniowa skala częstotliwości
w_vec = 2 * pi * f_vec; % Konwersja do pulsacji [rad/s]

% Obliczanie odpowiedzi częstotliwościowych
[mag_butt, ~] = freqresp(sys_butt, w_vec);
[mag_cheb1, ~] = freqresp(sys_cheb1, w_vec);
[mag_cheb2, ~] = freqresp(sys_cheb2, w_vec);
[mag_ellip, ~] = freqresp(sys_ellip, w_vec);

% Konwersja magnitudy do dB
mag_butt_dB = 20 * log10(abs(squeeze(mag_butt)));
mag_cheb1_dB = 20 * log10(abs(squeeze(mag_cheb1)));
mag_cheb2_dB = 20 * log10(abs(squeeze(mag_cheb2)));
mag_ellip_dB = 20 * log10(abs(squeeze(mag_ellip)));

% Rysowanie
figure;
plot(f_vec / 1000, mag_butt_dB, 'LineWidth', 1.5); hold on;
plot(f_vec / 1000, mag_cheb1_dB, 'LineWidth', 1.5);
plot(f_vec / 1000, mag_cheb2_dB, 'LineWidth', 1.5);
plot(f_vec / 1000, mag_ellip_dB, 'LineWidth', 1.5);

% Linie pomocnicze dla specyfikacji
plot([0 fp/1000], [-Rp -Rp], 'k--', 'LineWidth', 1); % Górna granica pasma przepustowego
plot([fp/1000 fp/1000], [-Rp -150], 'k--', 'LineWidth', 1); % Pionowa linia dla fp
plot([fsb/1000 fsb/1000], [-Rs -150], 'k--', 'LineWidth', 1); % Pionowa linia dla fsb
plot([fsb/1000 fs/1000], [-Rs -Rs], 'k--', 'LineWidth', 1); % Dolna granica pasma zaporowego

% Ustawienia wykresu
title('Charakterystyki Amplitudowe Filtrów Antyaliasingowych');
xlabel('Częstotliwość [kHz]');
ylabel('Amplituda [dB]');
legend(sprintf('Butterworth (n=%d)', n_butt), ...
       sprintf('Czebyszew I (n=%d)', n_cheb1), ...
       sprintf('Czebyszew II (n=%d)', n_cheb2), ...
       sprintf('Eliptyczny (n=%d)', n_ellip), ...
       'Wymagania', 'Location', 'southwest');
grid on;
axis([0 fs/1000 -100 5]); % Ograniczenie osi dla lepszej czytelności
hold off;

% Wybór najkorzystniejszego filtra
orders = [n_butt, n_cheb1, n_cheb2, n_ellip];
filter_names = {'Butterworth', 'Czebyszew I', 'Czebyszew II', 'Eliptyczny'};
min_order = min(orders);
best_filter_indices = find(orders == min_order);

fprintf('Najkorzystniejszy filtr (najniższy rząd spełniający wymagania):\n');
for i = 1:length(best_filter_indices)
    fprintf('- %s (rząd n = %d)\n', filter_names{best_filter_indices(i)}, min_order);
end
disp('Filtr eliptyczny zazwyczaj oferuje najniższy rząd dla zadanych specyfikacji,')
disp('ale kosztem zafalowań zarówno w paśmie przepustowym, jak i zaporowym.')
disp('Wybór zależy od priorytetów projektowych (rząd vs charakterystyka zafalowań).');
