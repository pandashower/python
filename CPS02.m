%DLA DOciekliwych w 2 i 4 zad do zrobienia i w 5 co dokładnie out!




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
y = A * x; 

% Rekonstrukcja (IDCT)
xs = S * y;  

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
mse_error = mean((x - xr).^2);
fprintf('Średni błąd kwadratowy (MSE) rekonstrukcji sygnału = %g\n', mse_error);


% zmiana f2 z 100 na 107
f2 = 107;
x2 = A1 * sin(2*pi*f1*t) + A2 * sin(2*pi*f2*t) + A3 * sin(2*pi*f3*t);
x2 = x2';
y2 = A * x2;

figure;
plot(f, y2);
grid on;

xrr = S * y2;

mse_error = mean((x - xrr).^2);
fprintf('Średni błąd kwadratowy (MSE) rekonstrukcji sygnału przy f2=107Hz = %g\n', mse_error);

f1 = f1 + 2.5;
f2 = f2 + 2.5;
f3 = f3 + 2.5;
x3 = A1 * sin(2*pi*f1*t) + A2 * sin(2*pi*f2*t) + A3 * sin(2*pi*f3*t);
x3 = x3.';
y3 = A * x3;

xrrr = S * y3;

figure;
plot(f, y3);
grid on;

mse_error = mean((x - xrrr).^2);
fprintf('Średni błąd kwadratowy (MSE) rekonstrukcji sygnału przy +2,5Hz = %g\n', mse_error);

% Tworzenie wektora indeksów k dla funkcji bazowych DCT
k = 10:10:30;
%k = 0:N-1; % żeby zaobserwować ze to dla wielokrotności piątki

% Obliczenie częstotliwości bazowych według wzoru:
f_k = (k / (2 * N)) * fs;
disp(" ");
disp('Nasze częstotliwości bazowe [Hz], które pokrywają się naturalnie w "y",');
disp('ponieważ jest on przekształceniem sygnału w zestaw współczynników,');
disp('które opisują, jak bardzo każda częstotliwość bazowa przyczynia się ');
disp('do jego budowy:');
disp(f_k);
disp("Dla k:")
disp(k)

% jak widać jeśli zmienie jedną z częstotliwości, np. na 107 to nie będzie 
% ona wielokrotnością 5 Hz i nie pokryje się z żadną częstotliwością bazową DCT.
% W takim przypadku DCT nie może jej dokładnie przedstawić za pomocą jednego
% współczynnika – energia tej częstotliwości "rozlewa się" na wiele sąsiednich
% współczynników. To zjawisko nazywa się przeciekiem częstotliwościowym 
% (ang. frequency leakage). W efekcie odwrotna transformacja (IDCT) nie 
% odtworzy sygnału idealnie, a błąd rekonstrukcji wzrośnie. 





%% Zad.5 (1,25p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

[x, fs] = audioread('mowa.wav'); % fs=8000Hz
plot(x); title('Oryginalny sygnał mowy'); 
xlabel('Próbka'); ylabel('Amplituda');
soundsc(x, fs); 
pause;

c = dct(x);

figure;
stem(c);
title('Współczynniki DCT');
xlabel('Indeks'); ylabel('Wartość');
pause;


% rekonstrukcja sygnału na bazie 25% pierwszych współczynników
N = length(c);
p25 = floor(0.25 * N);
c25 = [c(1 : p25); zeros(N - p25, 1)];  % wyzeruj pozostałe współczynniki
y25 = idct(c25);

figure;
plot(y25);
title('Rekonstrukcja z 25% współczynników DCT (pierwszych)');
xlabel('Próbka'); ylabel('Amplituda');
soundsc(y25, fs);
pause;

% rekonstrukcja na bazie 75% współczynników (np. ostatnich)
p75 = floor(0.75 * N);
c75 = [zeros(N - p75, 1); c(N - p75 + 1 : end)];
y75 = idct(c75);

figure;
plot(y75);
title('Rekonstrukcja z 75% współczynników DCT (ostatnich)');
xlabel('Próbka'); ylabel('Amplituda');
soundsc(y75, fs);
pause;



% (+0.25 pkt) Dodanie zakłócenia sinusoidalnego o częstotliwości 250 Hz
% UWAGA na orientację wektora (x musi być kolumną lub wierszem konsekwentnie).

% Dodaj sygnał sinusoidalny do oryginalnego sygnału mowy.
x_noisy = x + 0.5 * sin(2 * pi * 250 / fs * (0:length(x)-1)');

% zakłócony sygnał i go odsłuc
figure;
plot(x_noisy);
title('Sygnał mowy z zakłóceniem 250 Hz');
xlabel('Próbka'); ylabel('Amplituda');
soundsc(x_noisy, fs);
pause;

% Oblicz DCT sygnału zakłóconego.
c_noisy = dct(x_noisy);
figure;
stem(c_noisy);
title('Współczynniki DCT sygnału z zakłóceniem');
xlabel('Indeks'); ylabel('Wartość');
pause;

% Przykładowa metoda usuwania zakłócenia:
% Możemy wyzerować pewien wąski zakres indeksów "w okolicach" harmonicznego zaburzenia.
% (Najprostsze „ręczne” wyzerowanie jakiegoś zakresu, np. 100:200 – w praktyce
% należałoby ustalić dokładny indeks częstotliwości 250 Hz w DCT.)

c_noisy(100:200) = 0;

% Odwrócona transformacja DCT po "usunięciu" zakłócenia
y_clean = idct(c_noisy);

figure;
plot(y_clean);
title('Sygnał po usunięciu zakłócenia');
xlabel('Próbka'); ylabel('Amplituda');
soundsc(y_clean, fs);
pause;






%{
% 4. Dodanie zakłócenia sinusoidalnego 250 Hz
x_noise = x + 0.5*sin(2*pi*250/fs*(0:length(x)-1)');
figure; plot(x_noise); title('Sygnał z zakłóceniem');
% soundsc(x_noise, fs);


c_noise = dct(x_noise);
stem(c_noise)
c_noise(200:250) = 0;
x_clean = idct(c_noise);

% 6. Odsłuch i wyświetlenie wyników
figure; plot(x_clean); title('Oczyszczony sygnał');
% soundsc(x_clean, fs);
%}




% widoczny kształt (duże wartości na początku wykresu i „zgaśnięcie” 
% w dalszej części) wynika z tego, że sygnał mowy ma najwięcej 
% energii w niższych częstotliwościach.

% Duże skoki na początku mogą odpowiadać składowej stałej 
% (tzw. DC, pierwszy współczynnik) oraz wolnym (niskim) składowym
% częstotliwościowym, które zwykle zawierają najwięcej „energii” sygnału.
% W dalszej części współczynniki zazwyczaj mają mniejszą amplitudę i 
% odpowiadają wyższym częstotliwościom. Dlaczego „pierwsze 25%” 
% współczynników często brzmi „lepiej” niż „ostatnie 75%”? DCT „zbliża się” 
% do analizy częstotliwości – pierwsze współczynniki opisują głównie wolne
% zmiany (czyli składowe niskoczęstotliwościowe), a dopiero dalsze współczynniki
% reprezentują szybsze zmiany (wysokie częstotliwości). W mowie (i w większości
% sygnałów audio) większość energii jest przenoszona w niższych częstotliwościach. 
% Dlatego odcięcie początku (pierwszych współczynników) mocno zniekształca
% sygnał – bo tracimy bazową strukturę dźwięku.
