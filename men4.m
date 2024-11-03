clear all; close all; clc;

%% (*)PROBLEM 4.1
clear all; close all; clc;
% Czy macierz otrzymana w wyniku mnozenia Ahat=L*U jest jak oryginalna?
% czy funkcja z wydruku 4.1 daje poprawne wyniki dla dużych macierzy
% kwadratowych? TU:macierz 10x10 - nie jest taka sama jak org,
% dla 2 jeszcze tak zwykle

function [L,U] = myLU_orginalna(A)
[N, N] = size(A);
if(0) % prosciej, wolniej -------------------------------
    L = eye(N); U = zeros(N, N);
    for i = 1:N
        for j = i:N
            U(i, j) = A(i, j) - L(i, 1:i-1) * U(1:i-1, j);
        end
        for j = i+1:N
            L(j, i) = 1 / U(i, i) * (A(j, i) - L(j, 1:i-1) * U(1:i-1, i));
        end
    end
else % trudniej, szybciej -------------------------------
    U = A; L = eye(N);
    for i = 1:N-1
        for j = i+1:N
            L(j, i) = U(j, i) / U(i, i);
            U(j, i:N) = U(j, i:N) - L(j, i) * U(i, i:N);
        end
    end
end
end

% Testowanie dekompozycji LU
A = rand(10);  % Przykład z losową macierzą 10x10
[L, U] = myLU_orginalna(A);

% Sprawdzenie, czy A jest w przybliżeniu równe L*U
A_hat = L * U;
if norm(A - A_hat) <= 0
    disp('Dekompozycja LU jest poprawna. A jest równe L*U.');
else
    disp('Dekompozycja LU nie jest poprawna.');
end

%% Dodadaj do programu mozliwość obliczenia według wzoru (4.21). 
clear all; close all; clc;
% Porównaj z funkcja Matlaba [L,U]=lu(A). TU: macierz 10x10

function [L, U] = myLU(A)
    % Pobranie rozmiaru macierzy A
    [N, N] = size(A);
    
    % Inicjalizacja L i U zgodnie z (4.21)
    L = eye(N); 
    U = zeros(N, N);
    
    % Obliczanie elementów L i U przy użyciu formuł z równania (4.21)
    for i = 1:N
        for j = i:N
            % sumuje skalary zamiast wektorów
            U(i, j) = A(i, j) - sum(L(i, 1:i-1) .* U(1:i-1, j)');
        end
        for j = i+1:N
            % sumuje skalary zamiast wektorów
            L(j, i) = (A(j, i) - sum(L(j, 1:i-1) .* U(1:i-1, i)')) / U(i, i);
        end
    end
end

% Testowanie dekompozycji LU
A = rand(10);  % Przykład z losową macierzą 
[L, U] = myLU(A);

% Sprawdzenie, czy A jest w przybliżeniu równe L*U
A_hat = L * U;
if norm(A - A_hat) < 1e-10
    disp('Dekompozycja LU jest poprawna. A jest w przybliżeniu równa L*U.');
else
    disp('Dekompozycja LU nie jest poprawna.');
end

% Porównanie z wbudowaną funkcją lu w MATLAB
[L_builtin, U_builtin] = lu(A);
A_builtin_hat = L_builtin * U_builtin;

% Wyświetlenie różnic w celu porównania
disp('Różnica między moją dekompozycją LU a wbudowaną funkcją LU w MATLAB:');
disp(norm(A_hat - A_builtin_hat));





%% (*)PROBLEM 4.2
clear all; close all; clc;

function L = myCholesky(A)
    
    % Pobierz rozmiar macierzy A
    [N, ~] = size(A);
    
    % Inicjalizuj macierz L jako zerową macierz N x N
    L = zeros(N, N);
    
    % Obliczanie elementów L zgodnie ze wzorami
    for i = 1:N
        % Elementy przekątne
        L(i, i) = sqrt(A(i, i) - sum(L(i, 1:i-1).^2));
        
        % Elementy poniżej przekątnej
        for j = i+1:N
            L(j, i) = (A(j, i) - sum(L(j, 1:i-1) .* L(i, 1:i-1))) / L(i, i);
        end
    end
end

% Testowanie dekompozycji Choleskiego
A = rand(5); 
A = A' * A;  % Upewnienie się, że A jest symetryczna dodatnio określona

L = myCholesky(A);

% Sprawdzenie, czy A jest w przybliżeniu równe L * L'
A_hat = L * L';
if norm(A - A_hat) < 1e-10
    disp('Dekompozycja Choleskiego jest poprawna. A jest w przybliżeniu równa L * L''.');
else
    disp('Dekompozycja Choleskiego nie jest poprawna.');
end

% Porównanie z wbudowaną funkcją chol w MATLAB
L_builtin = chol(A, 'lower'); % Użycie opcji 'lower' dla dolnotrójkątnej wersji L

% Wyświetlenie różnic w celu porównania
disp('Różnica między moją dekompozycją Choleskiego a wbudowaną funkcją chol w MATLAB:');
disp(norm(L - L_builtin));





%% (*)(L)PROBLEM 4.9
clear all; close all; clc;
% L jest macierzą dolnotrójkątną, ale zawiera również jedynki na przekątnej.
function x = solveLU(A, b)
    % Użycie funkcji LU z macierzą permutacji
    [L, U, P] = lu(A);
    %A_back = P'*L*U

    % Przekształcenie prawej strony za pomocą permutacji
    b_permuted = P * b;
    
    % Liczba równań (rozmiar macierzy A)
    N = length(b);
    
    % Krok 1: Rozwiąż układ L * y = P * b (podstawianie w przód)
    y = zeros(N, 1);  % Inicjalizacja wektora y
    for i = 1:N
        y(i) = b_permuted(i) - L(i, 1:i-1) * y(1:i-1);
    end
    
    % Krok 2: Rozwiąż układ U * x = y (podstawianie wstecz)
    x = zeros(N, 1);  % Inicjalizacja wektora x
    for i = N:-1:1
        x(i) = (y(i) - U(i, i+1:N) * x(i+1:N)) / U(i, i);
    end
end

% Przykład użycia
A = [1, 2, 1; 3, 2, 1; 1, 1, 1];
b = [5; 17; 4];
x = solveLU(A, b);

disp('Rozwiązanie wektora x:');
disp(x);

% Sprawdzenie poprawności rozwiązania
b_calc = A * x;  % Obliczamy A * x
disp('Obliczony wektor b (A * x):');
disp(b_calc);

% Sprawdzenie różnicy między obliczonym b a oryginalnym b
difference = norm(b - b_calc);
disp('Różnica między obliczonym a oryginalnym b:');
disp(difference);

% Jeśli różnica jest bardzo mała, uznajemy rozwiązanie za poprawne
tolerance = 1e-10;  % Ustal tolerancję na małe błędy numeryczne
if difference < tolerance
    disp('Rozwiązanie jest poprawne.');
else
    disp('Rozwiązanie nie jest poprawne.');
end
%{
Macierz permutacji P wprowadza „przemieszczenie” wierszy, 
które pomaga poprawnie przeprowadzić dekompozycję LU, 
nawet jeśli w trakcie algorytmu pojawiłyby się zera na przekątnej,
co uniemożliwiłoby bezpośredni podział. Moje LU ze wzory zdaje się działać
w tym wypadku ale polecenie kazało użyć f. Matlaba
%}





%% (*)(*)(*)PROBLEM 4.5 - 1 cześć zadania
clear all; close all; clc;
% program (4.2)
A = [  1 2; ...
       3 4  ];
 
b = [    5; ...
        11  ];

x1 = inv(A)*b,          % x=A^(-1)*b
x2 = A\b;               % optymalne rozwiazywanie rown. Ax=b
%x3 = pinv(A)*b,        % x = inv( A'*A ) * A' * b , sprawdzisz?

bhat = A*x1;            % sprawdzenie
err = max(abs(x1-x2));  % blad

% Funkcja odwrot_rzad2: Obliczanie odwrotności macierzy 2x2 (4.31)
function [U] = odwrot_2x2(A)
    if det(A) == 0
        U = 0;
    else
        U = 1/(A(1,1)*A(2,2) - A(2,1)*A(1,2)) * [A(2,2) -A(1,2); -A(2,1) A(1,1)];
    end
end

% wywołanie odwrotu (4.31) i porównanie z wynikiem programu (4.2)
A_odwr = odwrot_2x2(A);
x1_prim = A_odwr*b,
if isequal(x1, x1_prim)
    disp("Odwrócenie prawidłowe");
end



%% 2 czesc zadania
clear all; close all; clc;
% Funkcja odwrot_rzad3: Obliczanie odwrotności macierzy 3x3 za pomocą minorów

clear all; close all; clc;
function [U] = odwrot_3x3(A)
    if det(A) == 0
        U = 0;
        disp("Wyznacznik nie może być równy 0");
    else
        N = size(A);
        Z = zeros(N(1), N(1));
        for i = 1:N(1)
            for j = 1:N(1)
                Matrix_dop = A;
                Matrix_dop(i,:) = [];      % Usunięcie i-tego wiersza
                Matrix_dop(:,j) = [];      % Usunięcie j-tej kolumny
                Z(i,j) = (-1)^(i+j) * det(Matrix_dop); % Obliczenie kofaktora
            end
        end
        U = (1 / det(A)) * Z.'; % Odwrócenie przy użyciu transpozycji macierzy dołączonej
    end
end

% program dla macierzy 3x3
A = [1 21 3;
     4 51 6;
     7  8 9];
b = [5; 11; 3];

A_odwr = odwrot_3x3(A),
A_z_inv = inv(A),

% Porównanie wyników
if isequal(round(A_odwr, 10), round(A_z_inv, 10)) % Porównanie z zaokrągleniem,
    % aby uwzględnić błędy numeryczne
    disp("Macierze odwrotne są równe w granicach błedu 10m. po , - odwrotność obliczona poprawnie");
else
    disp("Macierze odwrotne nie są równe - coś nie działa poprawnie");
end



%% 3 czesc zadania
clear all; close all; clc;
% Funkcja odwrot_recursive: Rekurencyjne obliczanie macierzy odwrotnej dowolnego rzędu

function U = odwrot_recursive(A, row_index, column_index, matrix_size, Z)
    if row_index < 1 || column_index < 1
        U = (1 / det(A)) * Z.'; % Zwróć wynik tylko, gdy macierz dopełnienia jest w pełni obliczona
        return;
    end

    Matrix_dop = create_dop_matrix(A, row_index, column_index);
    Z(row_index, column_index) = (-1)^(row_index + column_index) * det(Matrix_dop);

    if column_index > 1
        U = odwrot_recursive(A, row_index, column_index - 1, matrix_size, Z);
    elseif row_index > 1
        U = odwrot_recursive(A, row_index - 1, matrix_size, matrix_size, Z);
    else
        U = (1 / det(A)) * Z.'; % Zwróć wynik po zakończeniu rekurencji
    end
end

% Funkcja create_dop_matrix: Tworzy macierz dopełnienia
function Dop = create_dop_matrix(A, row_index, column_index)
    Dop = A;
    Dop(row_index, :) = []; % Usunięcie i-tego wiersza
    Dop(:, column_index) = []; % Usunięcie j-tej kolumny
end

% Główna funkcja dla macierzy dowolnego rzędu
A = [1 2 3 4; 3 4 6 8; 9 10 2 12; 9 11 15 3];
N = size(A);
Z = zeros(N(1), N(1));
A_inv = odwrot_recursive(A, N(1), N(1), N(1), Z);

% Sprawdzenie poprawności odwrotności
Z_corr = zeros(N(1), N(1));
A_corr = odwrot_recursive(A_inv, N(1), N(1), N(1), Z_corr);

if isequal(round(A, 10), round(A_corr, 10))
    disp("Wszystko zadziałało poprawnie");
else
    disp("Odwrócenie nie jest prawidłowe");
end

