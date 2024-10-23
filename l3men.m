clear all; close all; clc;

%% (*)(*)PROBLEM 3.8
clear all; close all; clc;
%Odsłanianie kanarka - usuwanie ryku słonia

N = 2^13; %ilość próbek 8196

% Wczytanie plików audio
[elephant, fs_e] = audioread('Elephant.wav', [1,N]); % Plik dźwiękowy słonia
[canary, fs_c] = audioread('Canary.wav', [1,N]); % Plik dźwiękowy kanarka

% Wykres kanarka i słonia
figure; plot(canary); title('kanarek'); xlabel('Czas'); ylabel('Amplituda');
figure; plot(elephant); title('slon'); xlabel('Czas'); ylabel('Amplituda');

% Suma sygnałów
combined_signal = elephant + canary;

% Wyświetlenie zsumowanego sygnału
figure; plot(combined_signal); title('Zsumowany sygnał (słoń + kanarek)');
xlabel('Czas'); ylabel('Amplituda');

% Widma częstotliwościowe
n=0:N-1; k=0:N-1;
A = sqrt(2/N)*cos( pi/N *(k'*n)); %specjalna macierz ortagonalna DCT
canary_spectrum = A * canary; % y = Ax
elephant_spectrum = A * elephant;
%figure; plot(canary_spectrum); title('Widmo orginalnego kanarka');
%figure; plot(elephant_spectrum); title('elephant freq spec');

% Widmo połączonego sygnału
combined_signal_fs = A * combined_signal;
figure; plot(combined_signal_fs); title('Widmo combined signal (słoń+kanarek)');
% odejmuje od laczonego widma, widmo kanarka
odsl_canary_fs = combined_signal_fs - elephant_spectrum;
%figure; plot(back_canary_fs); title('Widmo odsłoniętego kanarka');

% Transformacja odwrotna - SYNTEZA
odsl_canary    = A' * odsl_canary_fs; % dla DCT odwrotna transformacja to macierz
% transponowana, dla macierzy ortagonalnej o wart. R T^-1=T'
% częściowo tracimy te częstotliwości, które się dublowały ze słoniem 
figure; plot(odsl_canary); title('Odsłonięty kanarek');

% Odsłuchanie oryginalnego i oczyszczonego sygnału
sound(odsl_canary, fs_e); % dzwiek odsl kanarka
pause(length(odsl_canary)/fs_e + 1); % Pauza, aby poczekać na zakończenie odtwarzania



%% (*)(*)(*)PROBLEM 3.10
clear all; close all; clc;
load('babia_gora.dat'); % Wczytywanie danych
X = babia_gora; % Zakładam, że dane z pliku są w tej samej formie
size(X);

% 2D Rotation Section 
% Zmienne P = [x, y] dla obrotu w 2D
P_2D = X(:, 1:2); % Wybieramy pierwsze dwie kolumny (x, y) dla obrotu w 2D

% Definiowanie kąta obrotu dla 2D (obrót wokół osi Z)
theta = 45 / 180 * pi; % Kąt 45 stopni (można zmieniać)

% Macierz obrotu w 2D
R_2D = [cos(theta), -sin(theta);
        sin(theta), cos(theta)];

% Obrót punktów 2D
P_rot_2D = (R_2D * P_2D')'; % Transponowanie po obrocie, aby dopasować do formatu

% Wizualizacja przed i po obrocie w 2D
figure;
subplot(1, 2, 1); % Dzielimy figurę na dwie części (przed i po obrocie)
plot(P_2D(:, 1), P_2D(:, 2), 'b.');
title('Punkty przed obrotem 2D');
xlabel('X'); ylabel('Y');
grid on;

subplot(1, 2, 2);
plot(P_rot_2D(:, 1), P_rot_2D(:, 2), 'r.');
title('Punkty po obrocie 2D');
xlabel('X'); ylabel('Y');
grid on;

%%%%%%%%%%%%%%%
% Rysowanie punktów przed interpolacją
figure;
plot3(X(:,1), X(:,2), X(:,3), 'b.');
title('Wizualizacja plot3 przed obrotem');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;


% 3D Rotation Section
% Interpolacja danych
x = X(:,1); y = X(:,2); z = X(:,3);
vmin = min( min(x), min(y) ); % min zakres
vmax = max( max(x), max(y) ); % max zakres
[xi,yi] = meshgrid( vmin : (vmax-vmin)/200 : vmax ); % dopasowanie zakresu siatki
zi = griddata( x, y, z, xi, yi, 'linear' ); % interpolacja


% Definiowanie macierzy rotacji
ax = 180/180*pi; % Obrót wokół osi X
ay = 0/180*pi; % Obrót wokół osi Y
az = 0/180*pi; % Obrót wokół osi Z

Rx = [ 1, 0, 0;
       0, cos(ax), -sin(ax);
       0, sin(ax), cos(ax) ];% Macierz rotacji wokół osi X

Ry = [ cos(ay), 0, -sin(ay);
       0, 1, 0;
       sin(ay), 0, cos(ay) ]; % Macierz rotacji wokół osi Y

Rz = [ cos(az), -sin(az), 0;
       sin(az), cos(az), 0;
       0, 0, 1 ];           % Macierz rotacji wokół osi Z

% Przeprowadzanie obrotu (3 rotacje po kolei)
XR = Rz * Ry * Rx * X'; % Transponowanie X, aby pasowało do macierzy rotacji

% Rysowanie wyników po obrocie
figure;
plot3( XR(1,:), XR(2,:), XR(3,:), 'b.' );
title('Wizualizacja plot3 po obrocie');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;





%% (*)(*)PROBLEM 3.2
clear all; close all; clc;

R1 = 10; R2 = 20; R3 = 30; R0 = 40;
E1 = 1;  E2 = 2;  E3 = 3; E4 = 4;

A = [  R1+R2,   -R2,     0;  ...
         -R2, R2+R3,   -R3;  ...
           0,   -R3, R3+R0   ],
b = [ E1-E2; ...
      E2-E3; ...
        E3  ],

%Wersja 2
%{
A = [  R1+R2,   -R2,     0;  ...
         -R2, R2+R3,   -R3;  ...
           0,   -R3, R3+R0   ],
b = [ E1-E2; ...
      E2-E3; ...
        E3-E4  ],
%}
% x = ?

x1 = inv(A)*b;   % inv(A)  = A^(-1)
x2 = pinv(A)*b;  % pinv(A) = (A^T * A)^(-1) * A^T
x3 = A \ b;      % minimaliacja bledu sredniokwadratowego

% Metoda Cramera
for k=1:length(b)
    Ak = A; Ak(:,k) = b; % (w,k) = (:,k)
    x4(k) = det( Ak ) / det(A); 
end    
x4 = x4.';
[ x1, x2, x3, x4 ]
