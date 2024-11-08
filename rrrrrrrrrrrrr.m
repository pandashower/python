clear all; close all; clc;

% Parametry terenu i stacji referencyjnych
width = 80; % szerokość terenu (m)
height = 70; % wysokość terenu (m)
num_robots = 90; % liczba robotów

% Pozycje stacji referencyjnych
stations = [0, 0; 
            width, 0; 
            width, height; 
            0, height];

% Generowanie losowych pozycji robotów
robots_pos = [width * rand(num_robots, 1), height * rand(num_robots, 1)];

% Zmienna przechowująca oszacowane pozycje robotów
estimated_pos = zeros(num_robots, 2);

odchylenie_standardowe = 2; % odchylenie standardowe szumu (deg)

% Obliczanie oszacowanych pozycji robotów
for i = 1:num_robots
    A = [];
    b = [];
    for j = 1:size(stations, 1)
        % Prawdziwy kąt nadejścia sygnału
        true_angle = atand((robots_pos(i, 2) - stations(j, 2)) / (robots_pos(i, 1) - stations(j, 1)));
        
        % Szum
        angle_noise = odchylenie_standardowe * randn;
        
        % Oszacowany kąt nadejścia sygnału
        estimated_angle = true_angle + angle_noise;
        
        % Wyznaczanie równań dla algorytmu najmniejszych kwadratów
        A = [A; tand(estimated_angle), -1];
        b = [b; tand(estimated_angle) * stations(j, 1) - stations(j, 2)];
    end
    
    % Algorytm najmniejszych kwadratów z pseudoinwersją
    estimated_pos(i, :) = pinv(A' * A) * (A' * b);
end

% Obliczanie średniego błędu lokalizacji
location_errors = sqrt(sum((robots_pos - estimated_pos).^2, 2));
mean_location_error = mean(location_errors);

% Wyświetlanie wyników
figure;
hold on;
scatter(robots_pos(:, 1), robots_pos(:, 2), 'b', 'filled'); % prawdziwe pozycje robotów
scatter(estimated_pos(:, 1), estimated_pos(:, 2), 'r'); % oszacowane pozycje robotów
scatter(stations(:, 1), stations(:, 2), 'g', 'filled'); % pozycje stacji referencyjnych
legend('Prawdziwe pozycje', 'Oszacowane pozycje', 'Stacje referencyjne');
xlabel('X (m)');
ylabel('Y (m)');
title(['Średni błąd lokalizacji: ', num2str(mean_location_error), ' m']);
grid on;
axis equal;
hold off;

% -----------------------------------------
% Symulacja mapy błędu lokalizacji dla każdego punktu na siatce 1x1 m
% -----------------------------------------

% Rozdzielczość siatki (1 metr)
[x_grid, y_grid] = meshgrid(1:width, 1:height);

% Macierz do przechowywania średniego błędu lokalizacji
error_map = zeros(height, width);

% Liczba symulacji dla uśrednienia błędów
num_simulations = 50;

% Funkcja pomocnicza do obliczania błędu lokalizacji
calculate_error = @(true_pos, estimated_pos) sqrt(sum((true_pos - estimated_pos).^2));

% Przeprowadzanie symulacji dla każdej pozycji na siatce
for x = 1:width
    for y = 1:height
        errors = zeros(num_simulations, 1);
        
        for sim = 1:num_simulations
            true_pos = [x, y];
            A = [];
            b = [];
            
            for j = 1:size(stations, 1)
                % Prawdziwy kąt nadejścia sygnału
                true_angle = atand((true_pos(2) - stations(j, 2)) / (true_pos(1) - stations(j, 1)));
                
                % Szum
                angle_noise = odchylenie_standardowe * randn;
                
                % Oszacowany kąt nadejścia sygnału
                estimated_angle = true_angle + angle_noise;
                
                % Wyznaczanie równań dla algorytmu najmniejszych kwadratów
                A = [A; tand(estimated_angle), -1];
                b = [b; tand(estimated_angle) * stations(j, 1) - stations(j, 2)];
            end
            
            % Algorytm najmniejszych kwadratów z pseudoinwersją
            estimated_pos = pinv(A' * A) * (A' * b);
            
            % Obliczanie błędu lokalizacji
            errors(sim) = calculate_error(true_pos, estimated_pos');
        end
        
        % Średni błąd lokalizacji dla pozycji (x, y)
        error_map(y, x) = mean(errors);
    end
end

% Wizualizacja mapy błędów lokalizacji
figure;
pcolor(transpose(error_map));
shading interp;
colorbar;
title('Średni błąd lokalizacji w zależności od pozycji robota');
xlabel('X (m)');
ylabel('Y (m)');
axis equal tight;
