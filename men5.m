clear all; close all; clc;


%% (*)(*)PROBLEM 5.17
clear all; close all; clc;

% Wczytanie pliku z danymi topograficznymi
data = load('babia_gora.dat');
x = data(:,1); % współrzędne x
y = data(:,2); % współrzędne y
z = data(:,3); % wartości wysokości

figure; grid; plot3( data(:,1), data(:,2), data(:,3), 'b.' );


% Próbkowanie połowy danych do oceny dokładności interpolacji
half_x = data(1:2:end,1); half_y = data(1:2:end,2); half_z = data(1:2:end,3);
missing_data = data(2:2:end, :);


% Ustalenie siatki interpolacji
xvar = min(x) : (max(x)-min(x))/200 : max(x);
yvar = min(y) : (max(y)-min(y))/200 : max(y);
[Xi, Yi] = meshgrid(xvar, yvar); %siatka interpolacji xi, yi


% Interpolacja i wizualizacja dla różnych metod
methods = {'nearest', 'linear', 'cubic', 'v4'}; % 'v4' zamiast 'spline'

for i = 1:length(methods)
    figure;
    % Interpolacja za pomocą bieżącej metody
    out = griddata(x, y, z, Xi, Yi, methods{i});

    % Interpolacja na połowie danych dla oceny błędów
    half_out = griddata(half_x, half_y, half_z, Xi, Yi, methods{i});
    
    % Wyświetlenie interpolowanej powierzchni
    surf(Xi, Yi, out);
    title(['Metoda interpolacji: ', methods{i}]);
    xlabel('X');
    ylabel('Y');
    zlabel('Wysokość');



    % Obliczanie błędów
    RMSE = 0;
    for j = 1:size(missing_data, 1)
        % Znalezienie indeksów w interpolowanej siatce dla brakujących punktów
        [~, x_index] = min(abs(xvar - missing_data(j, 1)));
        [~, y_index] = min(abs(yvar - missing_data(j, 2)));
        
        % Wartości interpolowane i rzeczywiste
        interpolated_value = half_out(y_index, x_index);
        actual_value = out(y_index, x_index);

        % Obliczanie błędów
        abs_error = abs(interpolated_value - actual_value);
        RMSE = RMSE + abs_error^2;               % Błąd kwadratowy
    end
    
    % Średnia wartość błędów
    N = size(missing_data, 1);
    RMSE = sqrt(RMSE / N);

    % Wyświetlenie wyników błędu dla bieżącej metody
    disp(['Metoda: ', methods{i}]);
    disp(['Błąd średniokwadratowy (RMSE): ', num2str(RMSE)]);
end





%% (*)(*)(*)PROBLEM 5.15
clear all; close all; clc;

% Definiowanie danych wejściowych
x = [0, 1, 2, 3, 4];      % Węzły interpolacji
y = [0, 0.8415, 0.9093, 0.1411, -0.7568]; % Wartości funkcji w węzłach

n = length(x);

% Obliczanie odległości między węzłami
h = diff(x);

% Inicjalizacja macierzy i wektorów
A = zeros(4*(n-1));
b = zeros(4*(n-1), 1);

row = 1;

% Warunki interpolacji w węzłach
for i = 1:n-1
    % s_i(x_i) = y_i
    A(row, (i-1)*4+1) = 1;
    b(row) = y(i);
    row = row + 1;
    
    % s_i(x_{i+1}) = y_{i+1}
    dx = x(i+1) - x(i);
    A(row, (i-1)*4+1) = 1;
    A(row, (i-1)*4+2) = dx;
    A(row, (i-1)*4+3) = dx^2;
    A(row, (i-1)*4+4) = dx^3;
    b(row) = y(i+1);
    row = row + 1;
end


% Warunki ciągłości pierwszej i drugiej pochodnej w węzłach wewnętrznych
for i = 1:n-2
    % s_i'(x_{i+1}) = s_{i+1}'(x_{i+1})
    dx = x(i+1) - x(i);
    A(row, (i-1)*4+2) = 1;
    A(row, (i-1)*4+3) = 2*dx;
    A(row, (i-1)*4+4) = 3*dx^2;
    A(row, i*4+2) = -1;
    b(row) = 0;
    row = row + 1;
    
    % s_i''(x_{i+1}) = s_{i+1}''(x_{i+1})
    A(row, (i-1)*4+3) = 2;
    A(row, (i-1)*4+4) = 6*dx;
    A(row, i*4+3) = -2;
    b(row) = 0;
    row = row + 1;
end


% Warunki brzegowe (splajn naturalny)
% s_1''(x_1) = 0
A(row, 3) = 2;
b(row) = 0;
row = row + 1;

% s_{n-1}''(x_n) = 0
dx = x(n) - x(n-1);
A(row, (n-2)*4+3) = 2;
A(row, (n-2)*4+4) = 6*dx;
b(row) = 0;

% Rozwiązanie układu równań
coeffs = A\b;

% Wyodrębnienie współczynników a_i, b_i, c_i, d_i
a = coeffs(1:4:end);
b_coeffs = coeffs(2:4:end);
c = coeffs(3:4:end);
d = coeffs(4:4:end);

% Wyświetlenie współczynników
disp('Współczynniki a_i:');
disp(a');
disp('Współczynniki b_i:');
disp(b_coeffs');
disp('Współczynniki c_i:');
disp(c');
disp('Współczynniki d_i:');
disp(d');

% Sprawdzenie poprawności poprzez porównanie z funkcją spline MATLABa
xx = linspace(x(1), x(end), 100);
yy_manual = zeros(size(xx));

for i = 1:n-1
    idx = xx >= x(i) & xx <= x(i+1);
    dx = xx(idx) - x(i);
    yy_manual(idx) = a(i) + b_coeffs(i).*dx + c(i).*dx.^2 + d(i).*dx.^3;
end

% Wykorzystanie funkcji spline MATLABa
yy_builtin = spline(x, y, xx);

% Wykres porównawczy
figure;
plot(xx, yy_manual, 'r-', 'LineWidth', 2);
hold on;
plot(xx, yy_builtin, 'b--', 'LineWidth', 2);
plot(x, y, 'ko', 'MarkerFaceColor', 'k');
legend('Interpolacja manualna', 'Interpolacja spline MATLAB', 'Węzły');
title('Porównanie interpolacji kubicznej');
xlabel('x');
ylabel('y');
grid on;






%% (*)PROBLEM 5.3
% Przykład interpolacji okręgu na podstawie jego N punktów
clear all; close all;

N = 4; % stopień wielomianów od 6git
i = (0 : N)'; % zmienna "i" wielomianu w węzłach ("rzadka")
xi = cos(2 * pi / N * i); % wartości funkcji x = kosinus w węzłach
yi = sin(2 * pi / N * i); % wartości funkcji y = sinus w węzłach

% Sprawdzenie wartości
[i, xi, yi]

% Wygenerowanie i pokazanie macierzy Vandermonde'a
X = vander(i);

% Obliczenie wielomianów dla zmiennej x i y
ax = inv(X) * xi;
ay = inv(X) * yi;

% Dokładne wartości (gęstsze próbki)
id = (0 : 0.01 : N)';
xd = cos(2 * pi / N * id);
yd = sin(2 * pi / N * id);

% Interpolowane wartości za pomocą wielomianów
x_interp = polyval(ax, id);
y_interp = polyval(ay, id);

% Rysowanie wykresów
figure;
plot(xi, yi, 'ko', xd, yd, 'r-', x_interp, y_interp, 'b--');
xlabel('x'); ylabel('y'); title('y = f(x)'); axis square; grid on;
legend('Punkty węzłowe', 'Dokładny okrąg', 'Interpolowany okrąg');
pause;

figure;
plot(i, xi, 'ko', id, xd, 'r-', id, x_interp, 'b--');
xlabel('i'); ylabel('x'); title('x = f(i)'); grid on;
legend('Punkty węzłowe', 'Dokładne x', 'Interpolowane x');
pause;

% Obliczanie błędu między interpolacją a dokładnymi wartościami
error_x = abs(x_interp - xd);
error_y = abs(y_interp - yd);

% Średni błąd dla x i y
mean_error_x = mean(error_x);
mean_error_y = mean(error_y);

% Wyświetlanie wyników błędów
disp(['Średni błąd interpolacji dla x: ', num2str(mean_error_x)]);
disp(['Średni błąd interpolacji dla y: ', num2str(mean_error_y)]);

if mean_error_x < 1e-2
    disp("error x jest maly")
end
if mean_error_y < 1e-2
    disp ("error y jest mały")
end
