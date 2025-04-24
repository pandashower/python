%% Filtr cyfrowy IIR – realizacja zadania z polecenia
% Plik tworzy cyfrowy filtr IIR typu BP na podstawie prototypu
% Butterwortha zapisanego w butter.mat. Następnie bada jego
% charakterystyki, filtruje sygnał testowy i porównuje wyniki
% z wbudowaną funkcją filter().
clear; close all; clc;

%% 1. Parametry podstawowe
fs  = 16e3;             % częstotliwość próbkowania [Hz]
fc1 = 1189;             % dolna częstotliwość graniczna [Hz]
fc2 = 1229;             % górna częstotliwość graniczna [Hz]

%% 2. Wczytanie filtru analogowego
% plik butter.mat musi znajdować się w tej samej lokalizacji
% i zawierać zmienne: z – zera, p – bieguny, k – wzmocnienie
S = load('butter.mat');
[z, p, k] = deal(S.z, S.p, S.k);

% Równanie transmitancji H_a(s)
[ba, aa] = zp2tf(z, p, k);

%% 3. Konwersja biliniowa → filtr cyfrowy H(z) - powoduje warpowanie f
% bo 
[bd, ad] = bilinear(ba, aa, fs);

%% 4. Charakterystyka amplitudowo–częstotliwościowa
Nfft = 4096;                      % rozdzielczość FFT
f_vec = linspace(0, fs/2, Nfft); % oś częstotliwości w Hz
w_vec = 2*pi*f_vec;              % oś częstotliwości w rad/s

Ha = freqs(ba, aa, w_vec);       % analogowy
[Hd, fz] = freqz(bd, ad, Nfft, fs); % cyfrowy

figure('Name', 'Charakterystyki filtru');
plot(f_vec, 20*log10(abs(Ha)+eps), 'LineWidth',1.3); hold on;
plot(fz,   20*log10(abs(Hd)+eps), 'LineWidth',1.3);
xline(fc1, '--k', 'LineWidth',1);
xline(fc2, '--k', 'LineWidth',1);
grid on; xlabel('f [Hz]'); ylabel('|H| [dB]');
legend({'analogowy','cyfrowy','f_{c1}','f_{c2}'}, 'Location','Best');

title('Charakterystyka amplitudowa filtru BP – analog vs cyfrowy');

%% 5. Generacja sygnału testowego
T  = 1;                          % czas trwania [s]
t  = 0:1/fs:T-1/fs;              % wektor czasu
x  = sin(2*pi*1209*t) + sin(2*pi*1272*t);

%% 6. Filtracja – własna implementacja różniczkowa
Nb = numel(bd); Na = numel(ad);
y  = zeros(size(x));             % wyjście filtru

for n = 1:numel(x)
    % część odpowiadająca zerom (b)
    accB = 0;
    for k = 1:Nb
        if n-k+1 > 0
            accB = accB + bd(k) * x(n-k+1);
        end
    end
    % część odpowiadająca biegunom (a)
    accA = 0;
    for k = 2:Na
        if n-k+1 > 0
            accA = accA + ad(k) * y(n-k+1);
        end
    end
    % wyjście z równania różnicowego
    y(n) = (accB - accA) / ad(1);
end

%% 7. Filtracja wbudowaną funkcją filter() (kontrola)
y_ref = filter(bd, ad, x);

%% 8. Porównanie w dziedzinie czasu
figure('Name','Czas');
subplot(2,1,1);
plot(t, x);  grid on;
title('Sygnał oryginalny'); ylabel('Amplituda');
subplot(2,1,2);
plot(t, y); hold on; plot(t, y_ref, ':');
legend({'implementacja własna','filter()'});
 grid on; xlabel('t [s]');
title('Sygnał po filtracji');

%% 9. Porównanie w dziedzinie częstotliwości
X  = abs(fft(x, Nfft));
Y  = abs(fft(y, Nfft));
f_fft = (0:Nfft-1)*(fs/Nfft);

figure('Name','Widma');
plot(f_fft(1:Nfft/2), X(1:Nfft/2)/max(X), 'LineWidth',1.3); hold on;
plot(f_fft(1:Nfft/2), Y(1:Nfft/2)/max(Y), 'LineWidth',1.3);
legend({'przed filtracją','po filtracji'});
xlabel('f [Hz]'); ylabel('Znormalizowana amplituda'); grid on;
title('Widmo sygnału');

%% 10. Ocena zgodności obu metod filtracji
max_err = max(abs(y - y_ref));
fprintf('Maksymalny błąd między własną implementacją a filter(): %g\n', max_err);

%% 11. Uwagi:
%  • Różnice w częstotliwościach granicznych wynikają z transformacji
%    biliniowej (zjawisko warpowania częstotliwości). Aby dokładnie
%    odwzorować fc1 oraz fc2, należałoby zastosować prewarping.






%% *** OPCJA +0.25 pkt – PRE‑WARPING ***
T   = 1/fs;          % [s ] – okres próbkowania
figure('Name','Charakterystyki – bez i z pre‑warpingiem'); hold on;
plot(f_vec, 20*log10(abs(Ha)+eps),       'LineWidth',1.2); % H(s)
plot(f_vec, 20*log10(abs(Hd)+eps),       'LineWidth',1.2); % H(z)

% obliczenie zdokładnionych pulsacji analogowych - zgodnie ze wzorem
w1_pw = (2/T) * tan(pi*fc1/fs);    % ω₁₍w₎
w2_pw = (2/T) * tan(pi*fc2/fs);    % ω₂₍w₎

% Porząd filtru band‑pass (liczba biegunów)
N_bp = numel(p);            % rząd band‑pass (z butter.mat)
N_lp = N_bp/2;              % odpowiadający rząd prototypu low‑pass

% Projekt analogu z pre‑warpingiem (Butterworth BP, analog)
% Utrzymujemy ten sam rząd – różnią się jedynie częstotliwości graniczne.
[z_pw, p_pw, k_pw] = butter(N_lp, [w1_pw w2_pw], 'bandpass', 's');
[ba_pw, aa_pw]    = zp2tf(z_pw, p_pw, k_pw); % H_w(s)

% Konwersja biliniowa nowego prototypu (bez dodatkowego pre‑warpu – już jest)
[bd_pw, ad_pw] = bilinear(ba_pw, aa_pw, fs); % H_w(z)

% Widma analog/cyfrowy po pre‑warpingu
Ha_pw  = freqs(ba_pw, aa_pw, w_vec);
[Hd_pw, ~] = freqz(bd_pw, ad_pw, Nfft, fs);

% Dorysowanie charakterystyk pre‑warped
plot(f_vec, 20*log10(abs(Ha_pw)+eps), '--', 'LineWidth',1.2);
plot(f_vec, 20*log10(abs(Hd_pw)+eps), '--', 'LineWidth',1.2);

% Oznaczenia częstotliwości granicznych
xline(fc1, ':k', 'LineWidth',1);
xline(fc2, ':k', 'LineWidth',1);

xlabel('f [Hz]'); ylabel('|H| [dB]'); grid on;
legend({'H_a(s) – bez PW', ...
        'H_d(z) – bez PW', ...
        'H_{aw}(s) – pre‑warp', ...
        'H_{dw}(z) – pre‑warp', ...
        'f_{c1}','f_{c2}'}, 'Location','Best');
title({'Charakterystyka amplitudowa BP'; ...
       'porównanie: bez i z pre‑warpingiem'});

%% Komentarz:
% Pre‑warping przesuwa ω₁, ω₂ w analogu tak, by po BLT pasmo cyfrowe
% dokładnie pokrywało się z 1189 Hz i 1229 Hz.  Na wykresie widać, że
% czerwone przerywane H_{dw}(z) (z PW) przecina oś 0 dB niemal idealnie
% w miejscach linii f_{c1}, f_{c2}, podczas gdy filtr bez korekty (niebieski)
% jest delikatnie przesunięty w prawo.  Ostrzeżenia o "Matrix close to
% singular" powinny wyraźnie się zmniejszyć dzięki mniejszym wartościom
% współczynników po ponownym zaprojektowaniu filtru.
