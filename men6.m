clear all; close all; clc;

%% (*)PROBLEM 6.1
clear all; close all; clc;

x_data = [1; 2; 3; 4; 5];  % Wartości x
y_data = [2.2; 2.8; 3.6; 4.5; 5.1];  % Wartości y

% Tworzenie macierzy A oraz wektora b
A = [x_data, ones(size(x_data))];  % Macierz A: [x 1]
b = y_data;  % Wektor b

% Rozwiązanie za pomocą metody najmniejszych kwadratów
x = (A' * A) \ (A' * b);  % Wyznaczanie wektora x: [a; b]

% Wynik
a = x(1);  % Współczynnik nachylenia
b = x(2);  % Przesunięcie

% Wyświetlenie wyników
fprintf('Dopasowana prosta: y = %.2fx + %.2f\n', a, b);

% Rysowanie wyników
figure;
scatter(x_data, y_data, 'filled');  % Rysowanie punktów danych
hold on;
y_fit = A * x;  % Wyznaczanie wartości dopasowanej linii
plot(x_data, y_fit, 'r-', 'LineWidth', 2);  % Rysowanie dopasowanej linii
xlabel('x');
ylabel('y');
title('Dopasowanie linii prostej metodą najmniejszych kwadratów');
legend('Punkty danych', 'Dopasowana linia');
grid on;





%% (*)PROBLEM 6.7
clear all; close all; clc;

% definicja przedziału dla x
x = -1:0.01:1;

C = cell(1, 8); % komórki na wielomiany od C0 do C7
C{1} = ones(size(x)); % C0(x) = 1
C{2} = x; % C1(x) = x
C{3} = 2*x.*C{2}-C{1}; % C2 ...
C{4} = 2*x.*C{3}-C{2};
C{5} = 2*x.*C{4}-C{3};
C{6} = 2*x.*C{5}-C{4};
C{7} = 2*x.*C{6}-C{5};
C{8} = 2*x.*C{7}-C{6};

% wykres
figure;
hold on;
plot(x, C{1}, 'DisplayName', 'C_0(x)');
plot(x, C{2}, 'DisplayName', 'C_1(x)');
plot(x, C{3}, 'DisplayName', 'C_2(x)');
plot(x, C{4}, 'DisplayName', 'C_3(x)');
plot(x, C{5}, 'DisplayName', 'C_4(x)');
plot(x, C{6}, 'DisplayName', 'C_5(x)');
plot(x, C{7}, 'DisplayName', 'C_6(x)');
plot(x, C{8}, 'DisplayName', 'C_7(x)');

% legenda
legend show;
title('Wielomiany Czebyszewa C_0(x) do C_7(x) na przedziale [-1, 1]');
xlabel('x');
ylabel('C_n(x)');
grid on;
hold off;





%% (*)PROBLEM 6.8
clear all; close all; clc;

x = -3:0.1:3;

erf_matlab = erf(x);
erf_taylorem = (2 / sqrt(pi)) * (x - (x.^3) / 3 + (x.^5) / 10 - (x.^7) / 42);

% wykres
figure;
plot(x, erf_matlab, 'b-', 'LineWidth', 1.5); % wykres Matlab
hold on;
plot(x, erf_taylorem, 'r--', 'LineWidth', 1.5); % wykres Taylor

% ust osi i siatki
ylim([-1.2, 1.2]);
xlim([-2.5, 2.5]);
grid on;

title('Porównanie funkcji błędu (erf) i przybliżenia szeregiem Taylora', 'FontSize', 10);
xlabel('x', 'FontSize', 12);
ylabel('Wartość funkcji', 'FontSize', 12);
legend({'Funkcja erf Matlab', 'Przybliżenie szeregiem Taylora'}, 'Location', 'best', 'FontSize', 10);
hold off;


