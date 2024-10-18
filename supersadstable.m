clear all; close all; clc;
%jak c
% Parametry
Pt = 5; % Moc nadajnika w watach
f = 3e9; % Częstotliwość w Hz
c = 3e8; % Prędkość światła w m/s
lambda = c / f; % Długość fali
reflection_coeff = 0.7; % Współczynnik odbicia
x_bs = 110; % Pozycja stacji bazowej
y_bs = 190;
v = 30; % Prędkość użytkownika (m/s)

% Pozycje ścian
%{
walls = [
    60 40 200 40;  % Ściana 1
    70 90 180 90   % Ściana 2
];
%}
walls = [
    20 30 20 300;  % Ściana 1
    70 100 130 100   % Ściana 2
];

% Pozycja początkowa użytkownika
x_user = 30;
y_user_start = 10;

% Czas symulacji
t_total = 6; % sekundy
dt = 0.01; % krok czasowy (10 ms)
time = 0:dt:t_total; % Wektor czasu

% Trajektoria użytkownika
y_user = y_user_start + v * time;

% Odległość bezpośrednia (od stacji bazowej do użytkownika)
d_direct = sqrt((x_bs - x_user).^2 + (y_bs - y_user).^2);

% Obliczanie mocy sygnału bezpośredniego
Pr_direct = Pt * (lambda ./ (4 * pi * d_direct)).^2;

% Rysowanie wykresu okolicy z ścianami
figure;
hold on;
plot([walls(1,1), walls(1,3)], [walls(1,2), walls(1,4)], 'k-', 'LineWidth', 2); % Ściana 1
plot([walls(2,1), walls(2,3)], [walls(2,2), walls(2,4)], 'k-', 'LineWidth', 2); % Ściana 2
plot(x_bs, y_bs, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Stacja bazowa
plot(x_user * ones(size(y_user)), y_user, 'bo'); % Trajektoria użytkownika
title('Okolica ze ścianami i trajektorią użytkownika');
xlabel('X [m]');
ylabel('Y [m]');
grid on;

% Moc sygnału odbitego od ścian
Pr_reflected = zeros(size(time)); % Inicjalizacja

% Poprawiona pętla dla odbić od ścian
for i = 1:size(walls, 1)
    % Środek każdej ściany (upraszczamy odbicia od środka ściany)
    x_wall = (walls(i,1) + walls(i,3)) / 2;
    y_wall = (walls(i,2) + walls(i,4)) / 2;
    
    % Odległość stacja -> ściana -> użytkownik
    d_bs_wall = sqrt((x_bs - x_wall)^2 + (y_bs - y_wall)^2);
    d_wall_user = sqrt((x_wall - x_user).^2 + (y_wall - y_user).^2);
    
    % Całkowita odległość
    d_reflected = d_bs_wall + d_wall_user;
    
    % Moc odbita
    Pr_reflected = Pr_reflected + reflection_coeff * Pt * (lambda ./ (4 * pi * d_reflected)).^2;
end

% Sumaryczna moc odbierana
Pr_total = Pr_direct + Pr_reflected;

% Wykres mocy sygnału w funkcji czasu
figure;
plot(time, 10*log10(Pr_total)); % Moc w dBm
title('Moc sygnału odbieranego w funkcji czasu');
xlabel('Czas [s]');
ylabel('Moc [dBm]');
grid on;
