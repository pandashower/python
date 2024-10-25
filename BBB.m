clearvars;close all;clc;

% Stałe i dane
c = 3e8; % prędkość światła [m/s]
f = 6e9; % częstotliwość [Hz]
lambda = c / f; % długość fali [m]
Pt = 20e-3; % moc rutera [W]
Pt1 = Pt / 2; % moc w scenariuszu 1 (TDM)
Pt2 = Pt / 4; % moc w scenariuszu 2 (beamforming)
noise_power = 10^(-130/10); % moc szumów [W]

% Pozycje
router_pos = [30, 20]; % pozycja rutera
antenna1_pos = [29.9875, 20]; % pozycja anteny 1
antenna2_pos = [30.0125, 20]; % pozycja anteny 2
user1_pos = [100, 100]; % pozycja użytkownika 1
user2_pos = [140, 0]; % pozycja użytkownika 2

% Obliczenia dla ćwiczenia 1 (scenariusz 1: TDM)

% Odległości do użytkowników
d1_antenna1 = sqrt((user1_pos(1) - antenna1_pos(1))^2 + (user1_pos(2) - antenna1_pos(2))^2);
d1_antenna2 = sqrt((user1_pos(1) - antenna2_pos(1))^2 + (user1_pos(2) - antenna2_pos(2))^2);
d2_antenna1 = sqrt((user2_pos(1) - antenna1_pos(1))^2 + (user2_pos(2) - antenna1_pos(2))^2);
d2_antenna2 = sqrt((user2_pos(1) - antenna2_pos(1))^2 + (user2_pos(2) - antenna2_pos(2))^2);

% Tłumienie w przestrzeni swobodnej
L1_antenna1 = (lambda / (4 * pi * d1_antenna1))^2;
L1_antenna2 = (lambda / (4 * pi * d1_antenna2))^2;
L2_antenna1 = (lambda / (4 * pi * d2_antenna1))^2;
L2_antenna2 = (lambda / (4 * pi * d2_antenna2))^2;

% Sygnały odbierane w scenariuszu 1
Pr1 = Pt1 * (L1_antenna1 + L1_antenna2); % moc użytkownika 1
Pr2 = Pt1 * (L2_antenna1 + L2_antenna2); % moc użytkownika 2

% Obliczanie SNR
SNR1 = 10 * log10(Pr1 / noise_power);
SNR2 = 10 * log10(Pr2 / noise_power);

fprintf('SNR dla użytkownika 1 (scenariusz 1): %.2f dB\n', SNR1);
fprintf('SNR dla użytkownika 2 (scenariusz 1): %.2f dB\n', SNR2);

% Obliczenia dla ćwiczenia 2 (scenariusz 2: Orthogonal Beamforming)

% Przygotowanie faz (odwrotność faz)
Pr1_beamforming = Pt2 * (L1_antenna1 + L1_antenna2); % moc użytkownika 1
Pr2_beamforming = Pt2 * (L2_antenna1 + L2_antenna2); % moc użytkownika 2

% Obliczanie SNR w scenariuszu 2
SNR1_beamforming = 10 * log10(Pr1_beamforming / noise_power);
SNR2_beamforming = 10 * log10(Pr2_beamforming / noise_power);

fprintf('SNR dla użytkownika 1 (scenariusz 2): %.2f dB\n', SNR1_beamforming);
fprintf('SNR dla użytkownika 2 (scenariusz 2): %.2f dB\n', SNR2_beamforming);

% Ćwiczenie 3: użytkownicy w ruchu
v1 = 20; % prędkość użytkownika 1 [m/s]
v2 = 30; % prędkość użytkownika 2 [m/s]
t_end = 5; % czas symulacji [s]
dt = 0.1; % krok czasowy [s]

time = 0:dt:t_end; % przedział czasu
SNR1_moving = zeros(1, length(time));
SNR2_moving = zeros(1, length(time));

for t_idx = 1:length(time)
    t = time(t_idx);
    
    % Nowa pozycja użytkowników
    user1_new_pos = [100, 100 - v1 * t];
    user2_new_pos = [140, v2 * t];
    
    % Nowe odległości do anten
    d1_antenna1_new = sqrt((user1_new_pos(1) - antenna1_pos(1))^2 + (user1_new_pos(2) - antenna1_pos(2))^2);
    d1_antenna2_new = sqrt((user1_new_pos(1) - antenna2_pos(1))^2 + (user1_new_pos(2) - antenna2_pos(2))^2);
    d2_antenna1_new = sqrt((user2_new_pos(1) - antenna1_pos(1))^2 + (user2_new_pos(2) - antenna1_pos(2))^2);
    d2_antenna2_new = sqrt((user2_new_pos(1) - antenna2_pos(1))^2 + (user2_new_pos(2) - antenna2_pos(2))^2);
    
    % Tłumienie
    L1_antenna1_new = (lambda / (4 * pi * d1_antenna1_new))^2;
    L1_antenna2_new = (lambda / (4 * pi * d1_antenna2_new))^2;
    L2_antenna1_new = (lambda / (4 * pi * d2_antenna1_new))^2;
    L2_antenna2_new = (lambda / (4 * pi * d2_antenna2_new))^2;
    
    % Obliczanie mocy odbieranej
    Pr1_new = Pt2 * (L1_antenna1_new + L1_antenna2_new);
    Pr2_new = Pt2 * (L2_antenna1_new + L2_antenna2_new);
    
    % SNR
    SNR1_moving(t_idx) = 10 * log10(Pr1_new / noise_power);
    SNR2_moving(t_idx) = 10 * log10(Pr2_new / noise_power);
end

% Wykresy
figure;
plot(time, SNR1_moving, 'b', 'LineWidth', 2);
hold on;
plot(time, SNR2_moving, 'r', 'LineWidth', 2);
xlabel('Czas [s]');
ylabel('SNR [dB]');
legend('Użytkownik 1', 'Użytkownik 2');
title('SNR użytkowników w ruchu przez 5 sekund');
grid on;
