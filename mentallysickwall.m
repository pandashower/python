 close all; clear all; clc;
% Techniki Bezprzewodowe Lab 3
% Zaniki sygnału w komunikacji mobilnej

% Techniki Bezprzewodowe Lab 3
% Zaniki sygnału w komunikacji mobilnej

% Dane wejściowe
c = 3e8;                % Prędkość światła [m/s]
f = 2e9;                % Częstotliwość [Hz]
lambda = c / f;         % Długość fali [m]
Pt = 5;                 % Moc stacji bazowej [W]

t = 0:0.01:6;           % Czas od 0 do 6 s co 10 ms
N = length(t);

% Obliczenie prędkości, aby użytkownik przebył odległość od y=10 m do y=200 m w 6 s
v = (200 - 10) / 6;     % Prędkość użytkownika [m/s]

% Pozycje stacji bazowej i użytkownika
x_bs = 140; y_bs = 50;              % Pozycja BS
x_u = 30 * ones(1, N);              % Użytkownik porusza się wzdłuż x = 30 m
y_u = 10 + v * t;                   % Pozycja y użytkownika

% Ściany
% Pierwsza ściana od (60, 40) do (200, 40)
x_wall1_start = 60; x_wall1_end = 200; y_wall1 = 40;
% Druga ściana od (180, 70) do (180, 90)
x_wall2 = 180; y_wall2_start = 70; y_wall2_end = 90;

% Współczynnik odbicia
Gamma = 0.7;

% Inicjalizacja tablic
Pr = zeros(1, N);        % Moc odbierana
A0 = zeros(1, N);        % Amplituda fali bezpośredniej
A_ref1 = zeros(1, N);    % Amplituda odbicia od pierwszej ściany
A_ref2 = zeros(1, N);    % Amplituda odbicia od drugiej ściany

for i = 1:N
    % Odległość i faza fali bezpośredniej
    d0 = sqrt((x_bs - x_u(i))^2 + (y_bs - y_u(i))^2);
    phi0 = -2 * pi * d0 / lambda;
    A0(i) = (1 / d0) * exp(1j * phi0);
    
    % Odbicie od pierwszej ściany (y = 40)
    y_mirror1 = 2 * y_wall1 - y_bs;    % Pozycja odbicia BS względem ściany
    % Sprawdzenie, czy punkt odbicia znajduje się na ścianie
    x_i1 = x_u(i) + (x_bs - x_u(i)) * (y_wall1 - y_u(i)) / (y_bs - y_u(i));
    if x_i1 >= x_wall1_start && x_i1 <= x_wall1_end
        d_ref1 = sqrt((x_u(i) - x_bs)^2 + (y_u(i) - y_mirror1)^2);
        phi_ref1 = -2 * pi * d_ref1 / lambda;
        A_ref1(i) = (Gamma / d_ref1) * exp(1j * phi_ref1);
    else
        A_ref1(i) = 0;
    end
    
    % Odbicie od drugiej ściany (x = 180)
    x_mirror2 = 2 * x_wall2 - x_bs;    % Pozycja odbicia BS względem ściany
    % Sprawdzenie, czy punkt odbicia znajduje się na ścianie
    y_i2 = y_u(i) + (y_bs - y_u(i)) * (x_wall2 - x_u(i)) / (x_bs - x_u(i));
    if y_i2 >= y_wall2_start && y_i2 <= y_wall2_end
        d_ref2 = sqrt((x_u(i) - x_mirror2)^2 + (y_u(i) - y_bs)^2);
        phi_ref2 = -2 * pi * d_ref2 / lambda;
        A_ref2(i) = (Gamma / d_ref2) * exp(1j * phi_ref2);
    else
        A_ref2(i) = 0;
    end
    
    % Suma amplitud
    A_total = A0(i) + A_ref1(i) + A_ref2(i);
    % Moc odbierana (proporcjonalna do kwadratu modułu amplitudy)
    Pr(i) = (abs(A_total))^2;
end

% Wykres mocy odbieranej
figure;
plot(t, 10*log10(Pr));
xlabel('Czas [s]');
ylabel('Moc odbierana [dB]');
title('Moc sygnału odbieranego przez użytkownika w funkcji czasu');
grid on;