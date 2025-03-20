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

% sprawdzanie ortonormlanosci
for i = 1 : N-1
    for j = i+1 : N
        dotVal = dot(A(i,:), A(j,:));  % iloczyn skalarny wiersza i i j
        if abs(dotVal) < 10^-14
            fprintf('Wiersze %d i %d: ortonormalne (OK)\n', i, j);
        else
            fprintf('Wiersze %d i %d: NIE ortonormalne (iloczyn skalarny=%.5g)\n', i, j, dotVal);
        end
    end
end

% to nad pętlą tylko zeby sie wyświtliło niżej to tu dałem
fprintf('Maksymalna różnica implemetacji matlab od mojej: %e\n', max_diff)






%% Zad.2 (1,25p) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
xr = S * y;  

% Porównanie x i x_odtw
reconstruction_error = norm(x - xr);
disp('Błąd rekonstrukcji sygnału:');
disp(reconstruction_error)
% Jeśli transformacja jest idealna, błąd powinien być bardzo mały (w granicach
% błędów numerycznych).


%% (Dla dociekliwych) Sprawdzenie ortonormalności macierzy losowej i jej odwrotności
clear all; close all; clc;

N = 20;
A = randn(N);

for i = 1 : N-1
    for j = i+1 : N
        dotVal = dot(A(i,:), A(j,:));  % iloczyn skalarny wiersza i i j
        if abs(dotVal) < 10^-14
            fprintf('Wiersze %d i %d: ortonormalne (OK)\n', i, j);
        else
            fprintf('Wiersze %d i %d: NIE ortonormalne (iloczyn skalarny=%.5g)\n', i, j, dotVal);
        end
    end
end

% sprawdzamy "ortonormalność" wierszy, tzn. czy norma każdego wiersza ~ 1
rowNorms = zeros(N,1);
for i = 1:N
    rowNorms(i) = norm(A(i,:));
end
disp('Normy wierszy losowej macierzy A:');
disp(rowNorms);

S=inv(A); % moge tak ale cieżej z typowym w DCT A' bo aby ortonormlaność -> A^-1=A^T
% wierse kolumny muszą być parami prostopadłe i znormowane do 1,
% transpozycja bardziej optymalna obliczeniowo

% Sprawdzenie, czy S * A = I (macierz jednostkowa)
I_check = A * S;
disp('Maksymalne odchylenie randA (poza przekątną) od macierzy jednostkowej dla S * A = I:');
maxx = max(max(abs(I_check - eye(N)))); % 1 maks daje nawiększe odchylenia w danych kolumnach
disp(maxx); %blad duzy gdyby macierz nieodwracalna, kol/wiersze prawie liniowo zależne itp.

% Analiza i synteza dowolnego sygnału losowego (podobnie jak wyżej)
x = randn(N,1);
y = A * x;   % analiza
x_s = S * y; % synteza

% Czy x_s == x? (błąd rekonstrukcji)
reconstruction_error_rand = norm(x - x_s);
disp('Błąd rekonstrukcji sygnału przy macierzy randn() - owszem można dla każdego nieosobliwego A,');
disp('ale taka macierz nie ma dobrego rozkładu energi np. koncentracji niskich częstotiwości ani');
disp('jasnej interpretacji (losowe kolumny/wiersze niewiele wnoszą w kategoriach jakie częstotliwości lub wzorce');
disp('zawiera sygnał, wartość błedu:');
disp(reconstruction_error_rand);

N = 20; % wymiar macierzy
A_spoiled = zeros(N,N); % rezerwujemy miejsce na macierz DCT-II

% budowa macierzy DCT z błędem k+0.25
for k = 0 : N-1
    if k == 0
        s_k = sqrt(1/N);
    else
        s_k = sqrt(2/N);
    end
    
    for n = 0 : N-1
        A_spoiled(k+1, n+1) = s_k * cos( (pi*(k+0.25)/N) * (n + 0.5) );
    end
end

I_check = A_spoiled * A_spoiled';
disp('(~0 zeby ortogonalna) Maksymalne odchylenie A z k+0.25 (poza przekątną) od macierzy jednostkowej dla A^T*A = I:');
maxx = max(max(abs(I_check - eye(N))));
disp(maxx);


% Analiza i rekonstrukcja sygnału SZUMOWEGO
x_rand = randn(N,1);       % sygnał losowy
y_rand = A_spoiled * x_rand;  % "analiza" w zepsutej DCT

% Próba rekonstrukcji: w prawidłowej DCT odwracamy przez transpozycję,
% ale tu już macierz nie jest ściśle ortonormalna.
x_recon_rand = A_spoiled' * y_rand;  

error_rand = norm(x_rand - x_recon_rand);
disp(['Błąd rekonstrukcji sygnału losowego: ', num2str(error_rand)]);

% Analiza i rekonstrukcja sygnału HARMONICZNEGO
n = (0:N-1).'; 
freq = 2; % przykładowa częstotliwość dyskretna
x_sin = sin(2*pi*freq*n/N);  

y_sin = A_spoiled * x_sin;  
x_recon_sin = A_spoiled' * y_sin;

error_sin = norm(x_sin - x_recon_sin);
disp(['Błąd rekonstrukcji sygnału sinusoidalnego: ', num2str(error_sin)]);

% Komentarze końcowe
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

% budowa macierzy DCT (A) i IDCT (S)
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
    pause;
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



% (+0.25p) Dodanie zakłócenia sinusoidalnego o częstotliwości 250 Hz
% do oryginalnego sygnału mowy.
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
N = length(x_noisy);
f = 250;
k = round((2 * N * f) / fs); % wyprowadzajac wzór na obliczanie częstotliwości 
% bazowych dostajemy wzór na indeks w DCT
% biore kilka dziesiąt okalających również bo energia sie "rozlewa"

c_noisy(2348:2468) = 0; % po k liczym na dole to think about though

% Odwrócona transformacja DCT po "usunięciu" zakłócenia
y_clean = idct(c_noisy);

figure;
plot(y_clean);
title('Sygnał po usunięciu zakłócenia');
xlabel('Próbka'); ylabel('Amplituda');
soundsc(y_clean, fs);
pause;



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


