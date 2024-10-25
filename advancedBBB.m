% Symulacja Beamformingu dla AP komunikującego się z dwoma użytkownikami
% Autor: [Twoje Imię i Nazwisko]
% Data: [Aktualna Data]

clear all; close all; clc;

%% Sprawdzenie dostępności pakietów
if ~license('test', 'Phased_Array_System_Toolbox')
    error('Pakiet Phased Array System Toolbox jest wymagany do uruchomienia tego kodu.');
end

%% Stałe i Parametry
c = physconst('LightSpeed'); % Prędkość światła [m/s]
f = 6e9; % Częstotliwość [Hz]
lambda = c / f; % Długość fali [m]

% Definicja pojedynczej anteny izotropowej
antenna = phased.IsotropicAntennaElement('FrequencyRange',[1e9 10e9]);

% Pozycje rutera i anten
router_pos = [30; 20; 0]; % Pozycja rutera [m]
antenna_spacing = lambda / 2; % Odległość między antenami [m]

% Pozycje anten w układzie (globalne współrzędne)
antenna_positions = [29.9875, 30.0125; % x
                     20,       20;      % y
                     0,        0];      % z

% Definicja układu antenowego
array = phased.ConformalArray('Element', antenna, ...
                              'ElementPosition', antenna_positions);

% Pozycje użytkowników
user1_pos = [100; 100; 0]; % Użytkownik 1 [m]
user2_pos = [140; 200; 0]; % Użytkownik 2 [m]

% Moc całkowita rutera
Pt_total = 20e-3; % [W]

% Poziom szumów na wejściu odbiorników
Pn = 10^(-130/10); % [W]

%% Ćwiczenie 1: Scenariusz 1 (TDM)
% Moc rutera dzielimy na 2
Pt_per_user = Pt_total / 2; % [W]

% Użytkownik 1
% Obliczanie kierunku do użytkownika 1
vec_user1 = user1_pos - router_pos;
[az1, el1, d_user1] = cart2sph(vec_user1(1), vec_user1(2), vec_user1(3));
az1 = rad2deg(az1);
el1 = rad2deg(el1);

% Obliczanie wektora sterującego dla użytkownika 1
steervec = phased.SteeringVector('SensorArray', array, 'PropagationSpeed', c);
weights_user1 = steervec(f, [az1; el1]);

% Obliczanie zysku anteny w kierunku użytkownika 1
[~, gain_user1] = pattern(array, f, az1, el1, 'Weights', weights_user1);

% Obliczanie mocy odbieranej przez użytkownika 1
G_linear_user1 = 10^(gain_user1/10);
Pr_user1 = Pt_per_user * G_linear_user1 * (lambda / (4 * pi * d_user1))^2;
SNR_user1 = Pr_user1 / Pn;
SNR_dB_user1 = 10 * log10(SNR_user1);

% Użytkownik 2
vec_user2 = user2_pos - router_pos;
[az2, el2, d_user2] = cart2sph(vec_user2(1), vec_user2(2), vec_user2(3));
az2 = rad2deg(az2);
el2 = rad2deg(el2);

% Obliczanie wektora sterującego dla użytkownika 2
weights_user2 = steervec(f, [az2; el2]);

% Obliczanie zysku anteny w kierunku użytkownika 2
[~, gain_user2] = pattern(array, f, az2, el2, 'Weights', weights_user2);

% Obliczanie mocy odbieranej przez użytkownika 2
G_linear_user2 = 10^(gain_user2/10);
Pr_user2 = Pt_per_user * G_linear_user2 * (lambda / (4 * pi * d_user2))^2;
SNR_user2 = Pr_user2 / Pn;
SNR_dB_user2 = 10 * log10(SNR_user2);

% Wyświetlenie wyników
fprintf('Ćwiczenie 1:\n');
fprintf('SNR dla użytkownika 1: %.2f dB\n', SNR_dB_user1);
fprintf('SNR dla użytkownika 2: %.2f dB\n\n', SNR_dB_user2);

%% Ćwiczenie 2: Scenariusz 2 (Orthogonal Beamforming)
% Moc rutera dzielimy na 4
Pt_per_signal = Pt_total / 2; % Moc na sygnał
Pt_per_antenna = Pt_per_signal / 2; % Moc na antenę dla każdego sygnału

% Obliczanie wektorów kanałów
d_user1_antennas = vecnorm(antenna_positions - user1_pos, 2, 1).';
d_user2_antennas = vecnorm(antenna_positions - user2_pos, 2, 1).';

h1 = (lambda ./ (4 * pi * d_user1_antennas)) .* ...
     exp(-1j * 2 * pi * d_user1_antennas / lambda);
h2 = (lambda ./ (4 * pi * d_user2_antennas)) .* ...
     exp(-1j * 2 * pi * d_user2_antennas / lambda);

% Obliczanie matrycy kanałów
H = [h1.'; h2.']; % Rozmiar 2x2

% Projektowanie matrycy precodingu (Zero-Forcing)
W = H' / (H * H');

% Normalizacja mocy sygnałów
W = W ./ vecnorm(W);

% Transmitowany sygnał
s = [sqrt(Pt_per_signal); sqrt(Pt_per_signal)]; % Sygnały dla użytkowników
x = W * s;

% Obliczanie mocy transmitowanej przez każdą antenę (kontrola mocy)
Pt_transmitted = sum(abs(x).^2);
if Pt_transmitted > Pt_total
    x = x * sqrt(Pt_total / Pt_transmitted);
end

% Obliczanie mocy odbieranej przez użytkowników
y = H * x;
Pr_users = abs(y).^2;
SNR_users = Pr_users / Pn;
SNR_dB_users = 10 * log10(SNR_users);

% Wyświetlenie wyników
fprintf('Ćwiczenie 2:\n');
fprintf('SNR dla użytkownika 1: %.2f dB\n', SNR_dB_users(1));
fprintf('SNR dla użytkownika 2: %.2f dB\n\n', SNR_dB_users(2));

%% Ćwiczenie 3: Scenariusz 2 z ruchem użytkowników
% Parametry symulacji
t_total = 5; % Czas całkowity [s]
dt = 0.1; % Krok czasowy [s]
t = 0:dt:t_total; % Wektor czasu

% Ruch użytkownika 1
v1 = 20; % Prędkość [m/s]
user1_start = [100; 100; 0];
user1_end = [100; 0; 0];
user1_direction = (user1_end - user1_start) / norm(user1_end - user1_start);

% Ruch użytkownika 2
v2 = 40; % Prędkość [m/s]
user2_start = [140; 200; 0];
user2_end = [140; 300; 0];
user2_direction = (user2_end - user2_start) / norm(user2_end - user2_start);

% Inicjalizacja
SNR_user1_time = zeros(length(t), 1);
SNR_user2_time = zeros(length(t), 1);

% Symulacja w czasie
for idx = 1:length(t)
    % Aktualne pozycje użytkowników
    d1 = v1 * t(idx);
    if d1 <= norm(user1_end - user1_start)
        user1_pos = user1_start + user1_direction * d1;
    else
        user1_pos = user1_end;
    end

    d2 = v2 * t(idx);
    if d2 <= norm(user2_end - user2_start)
        user2_pos = user2_start + user2_direction * d2;
    else
        user2_pos = user2_end;
    end

    % Obliczanie wektorów kanałów
    d_user1_antennas = vecnorm(antenna_positions - user1_pos, 2, 1).';
    d_user2_antennas = vecnorm(antenna_positions - user2_pos, 2, 1).';

    h1 = (lambda ./ (4 * pi * d_user1_antennas)) .* ...
         exp(-1j * 2 * pi * d_user1_antennas / lambda);
    h2 = (lambda ./ (4 * pi * d_user2_antennas)) .* ...
         exp(-1j * 2 * pi * d_user2_antennas / lambda);

    % Obliczanie matrycy kanałów
    H = [h1.'; h2.']; % Rozmiar 2x2

    % Projektowanie matrycy precodingu (Zero-Forcing)
    W = H' / (H * H');

    % Normalizacja mocy sygnałów
    W = W ./ vecnorm(W);

    % Transmitowany sygnał
    x = W * s;

    % Kontrola mocy transmitowanej
    Pt_transmitted = sum(abs(x).^2);
    if Pt_transmitted > Pt_total
        x = x * sqrt(Pt_total / Pt_transmitted);
    end

    % Obliczanie mocy odbieranej przez użytkowników
    y = H * x;
    Pr_users = abs(y).^2;
    SNR_users = Pr_users / Pn;

    SNR_user1_time(idx) = SNR_users(1);
    SNR_user2_time(idx) = SNR_users(2);
end

% Konwersja SNR do dB
SNR_dB_user1_time = 10 * log10(SNR_user1_time);
SNR_dB_user2_time = 10 * log10(SNR_user2_time);

% Wykres SNR w czasie
figure;
plot(t, SNR_dB_user1_time, 'b-', 'LineWidth', 2);
hold on;
plot(t, SNR_dB_user2_time, 'r-', 'LineWidth', 2);
xlabel('Czas [s]');
ylabel('SNR [dB]');
title('SNR na wejściach odbiorników w czasie');
legend('Użytkownik 1', 'Użytkownik 2');
grid on;
