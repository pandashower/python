%% Zad.1 (1p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

N = 20; % wymiar macierzy
A = zeros(N,N); % rezerwujemy miejsce na macierz DCT-II

% --- Budowa macierzy analizującej ---
for k = 0 : N-1
    if k == 0
        s_k = sqrt(1/N); % standardowe współczynniki normalizujące zapewnijące
        %ortonormalność bazy, takiej normalizacji używa matlab według
        %wikipedii angielskiej
    else
        s_k = sqrt(2/N);
    end
    
    for n = 0 : N-1
        A(k+1, n+1) = s_k * cos( (pi*k/N) * (n + 0.5) );
    end
end

% Porównanie z implemengtacją matlaba
B = dctmtx(N);
diff = A - B;
max_diff = max(abs(diff(:)));

for i = 1 : N-1
    for j = i+1 : N
        dotVal = dot(A(i,:), A(j,:));  % iloczyn skalarny wiersza i i j
        if abs(dotVal) < 10^-14
            fprintf('Wiersze %d i %d: ortogonalne (OK)\n', i, j);
        else
            fprintf('Wiersze %d i %d: NIE ortogonalne (iloczyn skalarny=%.5g)\n', i, j, dotVal);
        end
    end
end

% to nad pętlą tylko zeby sie wyświtliło niżej to tu dałem
fprintf('Maksymalna różnica implemetacji matlab od mojej: %e\n', max_diff)






%% Zad.2 (1,25???p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

N = 20; % wymiar macierzy
A = zeros(N,N); % rezerwujemy miejsce na macierz DCT-II

for k = 0 : N-1
    if k == 0
        s_k = sqrt(1/N);
    else
        s_k = sqrt(2/N);
    end
    
    for n = 0 : N-1
        A(k+1, n+1) = s_k * cos( (pi*k/N) * (n + 0.5) );
    end
end

% DCT 2 z założenia (w zależności od implementacji w MATLAB-ie) jest ortonormalna 
% albo prawie ortonormalna z odpowiednimi współczynnikami skali. 
% W idealnym przypadku macierz odwrotna to po prostu transpozycja.
D = A; %USUNAC LATER DELETE
% Macierz odwrotna (IDCT, synteza). 
% Jeżeli D jest ortonormalna, to S = D' będzie macierzą odwrotną.
S = A';
I = eye(N);

% Sprawdzenie, czy S * A = I (macierz jednostkowa)
I_check = A * S;
disp('Maksymalne odchylenie (poza przekątną) od macierzy jednostkowej dla S * A = I:');
maxx = max(max(abs(I_check - eye(N)))); % 1 maks daje nawiększe odchylenia w danych kolumnach
disp(maxx);

% Sprawdzenie rekonstrukcji sygnału (perfekcyjna rekonstrukcja?)
x = randn(N,1); 
% Transformacja sygnału x do dziedziny DCT
X = A * x;      

% Rekonstrukcja (IDCT)
xs = S * X;  

% Porównanie x i x_odtw
reconstruction_error = norm(x - xs);
disp('Błąd rekonstrukcji sygnału:');
disp(reconstruction_error)
% Jeśli transformacja jest idealna, błąd powinien być bardzo mały (w granicach
% błędów numerycznych).


%% (Dla dociekliwych) Sprawdzenie ortonormalności macierzy losowej i jej odwrotności
% Teraz wygenerujmy kwadratową macierz A (np. 20x20) za pomocą randn().
N2 = 20;
A = randn(N2);

% Sprawdźmy "ortonormalność" w sensie: czy kolumny (bądź wiersze) tworzą bazę
% ortonormalną. Aby macierz była ortonormalna, A'*A = I (dla kolumn ortonormalnych).
% Dla losowej macierzy randn() prawdopodobnie nie będzie ortonormalna, 
% ale możemy to po prostu zobaczyć.
orth_check = A' * A;

disp('Macierz A''*A (jeśli A byłaby ortonormalna, to byłaby zbliżona do I):');
disp(orth_check);

% Możemy też sprawdzić odchylenie od I:
error_orth = norm(orth_check - eye(N2));
disp(['Błąd ||A''*A - I|| = ', num2str(error_orth)]);

%% 6. Wyznaczenie macierzy odwrotnej i sprawdzenie A * inv(A) = I
S_inv = inv(A);        % Macierz odwrotna do A
I_check2 = A * S_inv;  % Sprawdzenie

disp('Macierz A * inv(A) (powinna być zbliżona do I):');
disp(I_check2);

error_inv = norm(I_check2 - eye(N2));
disp(['Błąd ||A * inv(A) - I|| = ', num2str(error_inv)]);

%% 7. "Zepsucie" macierzy A i obserwacja wpływu na odwrotność
% Przykładowo, zmieńmy drastycznie jeden z elementów A.
A(1,1) = A(1,1) + 1e3;  % Duża zmiana w pierwszym elemencie

% Ponownie wyznaczmy odwrotność
S_inv_corrupted = inv(A);

% Sprawdźmy, co się dzieje z iloczynem A * S_inv_corrupted
I_check3 = A * S_inv_corrupted;

disp('Macierz A * inv(A) po "zepsuciu" (powinna mocno odbiegać od I, w zależności od modyfikacji):');
disp(I_check3);

error_inv_corrupted = norm(I_check3 - eye(N2));
disp(['Błąd ||A * inv(A) - I|| po modyfikacji A = ', num2str(error_inv_corrupted)]);

%% Komentarze końcowe
% - W części z DCT/IDCT pokazujemy, że jeśli macierz D jest ortonormalna,
%   to jej transpozycja D' jest odwrotnością, co prowadzi do idealnej rekonstrukcji sygnału.
% - W części "dla dociekliwych" generujemy losową macierz A, sprawdzamy, 
%   czy jest ortonormalna (zwykle nie będzie), a następnie wyznaczamy i testujemy jej odwrotność.
% - "Zepsucie" macierzy A (modyfikacja jednego elementu) często dramatycznie
%   zmienia odwrotność, co ilustruje wrażliwość obliczeń na perturbacje.







%% Zad.3 (3p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

N = 100;             % liczba próbek sygnału
fs = 1000;           % częstotliwość próbkowania [Hz]
t = (0:N-1)/fs;      % wektor czasu (od 0 do (N-1)/fp)
f1 = 50;
f2 = 100;
f3 = 150;
A1 = 50;
A2 = 100;
A3 = 150;

% sygnał x jest sumą trzech sinusów o powyższych parametrach.
x = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t) + A3*sin(2*pi*f3*t);
x=x';

% Budowa macierzy DCT (A) i IDCT (S)
A = zeros(N,N); 
for k = 0 : N-1
    if k == 0
        s_k = sqrt(1/N);
    else
        s_k = sqrt(2/N);
    end
    
    for n = 0 : N-1
        A(k+1, n+1) = s_k * cos( (pi*k/N) * (n + 0.5) );
    end
end

S = A';

% wyświetlanie wierszy macierzy A i kolumn macierzy S w pętli
for i = 1:N
    fprintf('Wiersz %d macierzy A:\n', i);
    disp(A(i,:));         
    
    fprintf('Kolumna %d macierzy S:\n', i);
    disp(S(:,i));         
end

figure;
for i =1:N
hold on;
    plot(A(i,:),"b-x")
    plot(S(:, i),"r")
    %pause;
end

y = A * x;
f = (0:N-1) * fs / N / 2;
n=1:N; %bez skalowania na czest. 
figure;
plot(f, y);
grid on;

xr = S * y;
reconstruction_error = max(abs(x - xr));
fprintf('Maksymalny błąd rekonstrukcji sygnału (DCT->IDCT) = %g\n', ...
         reconstruction_error);

% zmiana f2 z 100 na 107
f2 = 105;
x = A1 * sin(2*pi*f1*t) + A2 * sin(2*pi*f2*t) + A3 * sin(2*pi*f3*t);
x = x';
y = A * x;

figure;
plot(f, y);
grid on;

xr = S * y;
reconstruction_error = max(abs(x - xr));
fprintf('Maksymalny błąd rekonstrukcji sygnału przy f2=107Hz = %g\n', ...
         reconstruction_error);

f1 = f1 + 2.5;
f2 = f2 + 2.5;
f3 = f3 + 2.5;
x = A1 * sin(2*pi*f1*t) + A2 * sin(2*pi*f2*t) + A3 * sin(2*pi*f3*t);
x = x.';
y = A * x;

figure;
plot(f, y);
grid on;



%% Zad.5 (1,25???p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
