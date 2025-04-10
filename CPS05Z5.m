% --- Czyszczenie środowiska ---
clc;            % Czyści okno poleceń
clear;          % Usuwa zmienne z przestrzeni roboczej
close all;      % Zamyka wszystkie otwarte okna wykresów

% =========================================================================
% === 1. PROJEKT FILTRA TESTOWEGO (96 MHz +/- 1 MHz) =====================
% =========================================================================
fprintf('--- Projektowanie FILTRA TESTOWEGO ---\n');

% --- Parametry filtra TESTOWEGO ---
fc_test = 96e6;      % Częstotliwość środkowa (Hz)
bw_pass_test = 2e6;  % Szerokość pasma przepustowego (Hz) (±1 MHz)

% --- Definicja częstotliwości granicznych TESTOWEGO ---
fp1_test = fc_test - bw_pass_test/2; % 95 MHz
fp2_test = fc_test + bw_pass_test/2; % 97 MHz
fs1_test = fp1_test - 1e6; % 94 MHz
fs2_test = fp2_test + 1e6; % 98 MHz

% --- Wymagania dotyczące tłumienia ---
Ap_test = 3;     % Maksymalne tłumienie w paśmie przepustowym (dB)
Ast_test = 40;   % Minimalne tłumienie w paśmie zaporowym (dB)

% --- Konwersja częstotliwości do pulsacji kątowej (rad/s) ---
Wp_test = 2 * pi * [fp1_test, fp2_test];
Ws_test = 2 * pi * [fs1_test, fs2_test];

% --- Obliczenie minimalnego rzędu filtra TESTOWEGO ---
[n_test, ~] = ellipord(Wp_test, Ws_test, Ap_test, Ast_test, 's');
fprintf('Filtr TESTOWY - Obliczony minimalny rząd filtra: n = %d\n', n_test);

% --- Projektowanie filtra analogowego TESTOWEGO ---
[b_test, a_test] = ellip(n_test, Ap_test, Ast_test, Wp_test, 's');

% --- Analiza charakterystyki częstotliwościowej TESTOWEGO ---
w_plot_test = linspace(2*pi*90e6, 2*pi*102e6, 2000); % Rozszerzony zakres: 90-102 MHz
[h_test, w_resp_test] = freqs(b_test, a_test, w_plot_test);
f_resp_test = w_resp_test / (2 * pi);

% --- Wyświetlanie charakterystyki częstotliwościowej TESTOWEGO ---
figure;
plot(f_resp_test / 1e6, 20*log10(abs(h_test)));
grid on;
title('Charakterystyka częstotliwościowa - FILTR TESTOWY (96 MHz \pm 1 MHz)');
xlabel('Częstotliwość [MHz]');
ylabel('Amplituda [dB]');
ylim([-Ast_test-20, 5]);
xlim([90, 102]); % Rozszerzony zakres osi X

% --- Zaznaczenie punktów charakterystycznych TESTOWEGO ---
hold on;
line([fp1_test, fp1_test]/1e6, ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
line([fp2_test, fp2_test]/1e6, ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
line([fs1_test, fs1_test]/1e6, ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
line([fs2_test, fs2_test]/1e6, ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
line(xlim, [-Ap_test, -Ap_test], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 1);
line(xlim, [-Ast_test, -Ast_test], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
legend('Charakterystyka filtra', ...
    sprintf('Granica pasma przep. (%.0f MHz)', fp1_test/1e6), ...
    sprintf('Granica pasma przep. (%.0f MHz)', fp2_test/1e6), ...
    sprintf('Granica pasma zapor. (%.0f MHz)', fs1_test/1e6), ...
    sprintf('Granica pasma zapor. (%.0f MHz)', fs2_test/1e6), ...
    sprintf('Tłumienie w paśmie przep. (-%.0f dB)', Ap_test), ...
    sprintf('Tłumienie w paśmie zapor. (-%.0f dB)', Ast_test), ...
    'Location', 'southwest');
hold off;
fprintf('Wygenerowano wykres dla FILTRA TESTOWEGO.\n\n');

% =========================================================================
% === 2. PROJEKT FILTRA DOCELOWEGO (96 MHz +/- 100 kHz) ==================
% =========================================================================
fprintf('--- Projektowanie FILTRA DOCELOWEGO ---\n');

% --- Parametry filtra DOCELOWEGO ---
fc = 96e6;      % Częstotliwość środkowa (Hz)
bw_pass = 200e3;  % Szerokość pasma przepustowego (Hz) (±100 kHz)
bw_stop = 2e6;    % Odstęp do sąsiednich stacji (Hz) (±1 MHz)

% --- Definicja częstotliwości granicznych DOCELOWEGO ---
fp1 = fc - bw_pass/2; % 95.9 MHz
fp2 = fc + bw_pass/2; % 96.1 MHz
fs1 = fc - bw_stop/2; % 95 MHz
fs2 = fc + bw_stop/2; % 97 MHz

% --- Wymagania dotyczące tłumienia ---
Ap = 3;     % Maksymalne tłumienie w paśmie przepustowym (dB)
Ast = 40;   % Minimalne tłumienie w paśmie zaporowym (dB)

% --- Konwersja częstotliwości do pulsacji kątowej (rad/s) ---
Wp = 2 * pi * [fp1, fp2];
Ws = 2 * pi * [fs1, fs2];

% --- Obliczenie minimalnego rzędu filtra DOCELOWEGO ---
[n, ~] = ellipord(Wp, Ws, Ap, Ast, 's');
fprintf('Filtr DOCELOWY - Obliczony minimalny rząd filtra: n = %d\n', n);

% --- Projektowanie filtra analogowego DOCELOWEGO ---
[b, a] = ellip(n, Ap, Ast, Wp, 's');

% --- Analiza charakterystyki częstotliwościowej DOCELOWEGO ---
w_plot = linspace(2*pi*90e6, 2*pi*102e6, 2000); % Rozszerzony zakres: 90-102 MHz
[h, w_resp] = freqs(b, a, w_plot);
f_resp = w_resp / (2 * pi);

% --- Wyświetlanie charakterystyki częstotliwościowej DOCELOWEGO ---
figure;
plot(f_resp / 1e6, 20*log10(abs(h)));
grid on;
title('Charakterystyka częstotliwościowa - FILTR DOCELOWY (96 MHz \pm 100 kHz)');
xlabel('Częstotliwość [MHz]');
ylabel('Amplituda [dB]');
ylim([-Ast-20, 5]);
xlim([90, 102]); % Rozszerzony zakres osi X

% --- Zaznaczenie punktów charakterystycznych DOCELOWEGO ---
hold on;
line([fp1, fp1]/1e6, ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
line([fp2, fp2]/1e6, ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
line([fs1, fs1]/1e6, ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
line([fs2, fs2]/1e6, ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
line(xlim, [-Ap, -Ap], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 1);
line(xlim, [-Ast, -Ast], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
legend('Charakterystyka filtra', ...
    sprintf('Granica pasma przep. (%.1f MHz)', fp1/1e6), ...
    sprintf('Granica pasma przep. (%.1f MHz)', fp2/1e6), ...
    sprintf('Granica pasma zapor. (%.0f MHz)', fs1/1e6), ...
    sprintf('Granica pasma zapor. (%.0f MHz)', fs2/1e6), ...
    sprintf('Tłumienie w paśmie przep. (-%.0f dB)', Ap), ...
    sprintf('Tłumienie w paśmie zapor. (-%.0f dB)', Ast), ...
    'Location', 'southwest');
hold off;
fprintf('Wygenerowano wykres dla FILTRA DOCELOWEGO.\n');

% --- Uwagi ---
disp('--- Uwagi ---');
disp('Jeśli charakterystyka filtra jest niezadowalająca:');
disp('- Zwiększ ręcznie rząd filtra (zmienna ''n'').');
disp('- Złagódź wymagania (np. zmniejsz Ast lub zwiększ Ap).');
disp('- Wypróbuj inny typ filtra (np. ''cheby1ord'', ''cheby2ord'', ''butterord'').');
