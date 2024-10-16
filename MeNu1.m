%W sumie to 2 rozdział
close all; clear all; clc; %clc to clear command window
% x=fi( value, sign, noBitsAll, noBitsFract ), x.bin,
%   value - wartos´c liczby rzeczywistej do zakodowania, ´
%   sign - 0=bez znaku, 1=ze znakiem w kodzie U2,
%   noBitsAll - liczba wszystkich uzytych bitów, ˙
%   noBitsFract - liczba bitów tylko czesci ułamkowej

%2.13, 2.15 i jakieś losowe daje 6*(*)


%(*)PROBLEM 2.1
%{
x=fi( 67, 1, 8, 0 ), x.bin, 
a=fi( 141,0,8,0), a.bin,
b=fi( 115,0,8,0), b.bin,
c=fi( 115,1,8,0), c.bin,
d=fi(-115,1,8,0), d.bin,

dzien_urodzin_8_bezznaku = fi(25, 0, 8, 0), dzien_urodzin_8_bezznaku.bin
dzien_urodzin_8 = fi(25, 1, 8, 0), dzien_urodzin_8.bin
dzien_urodzin_16_bezznaku = fi(25, 0, 16, 0), dzien_urodzin_16_bezznaku.bin
dzien_urodzin_16 = fi(25, 1, 16, 0), dzien_urodzin_16.bin
%}



%(*)(L)PROBLEM 2.4
%{
num2bitstr( single( (1+1/4)*2^(-124) )) == num2bitstr( single( -5.877472*10^(-38) )),

c = 299792458; %[m/s]

%Prędkość swiatła zapisana w trybie single (32b) i double (64b)
bit_single_c = num2bitstr(single(c)),
bit_double_c = num2bitstr(double(c)),

%Odtworzenie wartości z ich reprezentacji bitowych
function [liczba_dec] = bitstr2num(liczba_bitstr)
%Konwersja z użyciem wzorów z rys. 2.6
S = liczba_bitstr(1);
    %wyznaczenie K, L, F
    if( length(liczba_bitstr)==64 )
        e = liczba_bitstr(2 : 2+(11-1)); %exponent wykładnik
        f = liczba_bitstr(2+11 : 2+11+(52-1)); %fractionc mantysa
        L = length(f); %długość L - mantysy
        F = 1; %bo 1 jest dodawane do sumy
        for i = 1:L
            F = F + str2double(f(i)) * 2^(-i);
        end

        K = length(e); % liczba bitów wykładnika
        E_raw = 0;
        for i = 1:K
            E_raw = E_raw + str2double(e(i)) * 2^(K-i);
        end
        E = E_raw - (2^(K-1) - 1);

        wartosc = (-1)^S*F*2^E;

        liczba_dec = wartosc;

    elseif( length(liczba_bitstr)==32 )
        e = liczba_bitstr(2 : 2+(8-1)); %exponent wykładnik
        f = liczba_bitstr(2+8 : 2+8+(23-1)); %fractionc mantysa
        L = length(f); %długość L - mantysy
        F = 1; %bo 1 jest dodawane do sumy
        for i = 1:L
            F = F + str2double(f(i)) * 2^(-i);
        end

        K = length(e); % liczba bitów wykładnika
        E_raw = 0;
        for i = 1:K
            E_raw = E_raw + str2double(e(i)) * 2^(K-i);
        end
        E = E_raw - (2^(K-1) - 1);

        wartosc = (-1)^S*F*2^E;

        liczba_dec = wartosc;
    else    
    end
end

a = bitstr2num(bit_single_c);
b = bitstr2num(bit_double_c);
fprintf('Wartość pr. św: %.10f\n', c);
fprintf('Wartość single: %.10f\n', a);
fprintf('Wartość double: %.10f\n', b);
%}



%(*)(L)PROBLEM 2.11
%{
 x = 10 ^(-50);
 y = 10 ^200;
 z = 10 ^300;
 a = (x*y)/z
 b = x*(y/z)
 c = (x/z)*y
%}



%(*)(*)(L)PROBLEM 2.13
%{
%b>>4ac
a = 100;
b = 1000;
c = -100;

delta = b^2 - 4*a*c;

if delta < 0
    error('Równanie nie ma rozwiązań rzeczywistych.');
end

x1 = (-b - sqrt(delta)) / (2*a);
x2 = (-b + sqrt(delta)) / (2*a);

fprintf('              x1 = %.20f\n', x1); %w double 15-17cyfr znaczących
fprintf('              x2 = %.20f\n', x2); %zmiennoprzecin./dziesiętnie, nowa linia do 20po prze

if(abs(x1) >= abs(x2)) %wartości bezwzględne
   x2_prim = c / (a * x1); %3 punkt przepisu kulinarnego
   fprintf('Dokładniejsze x2 = %.20f\n', x2_prim);
   
   diff = x2 - x2_prim;
   fprintf('         Różnica = %.20f\n', diff);
else
   x1_prim = c / (a * x2);
   fprintf('Dokładniejsze x1 = %.20f\n', x1_prim);

   diff = x1 - x1_prim;
   fprintf('Różnica = %.20f\n', diff);
end

disp("rożnica: "+diff)
%}



%(*)(*)(*)PROBLEM 2.15
%format long;

czestotliwosc = 40;  % Częstotliwość sygnału w Hz - liczba powt. sinsuoidy na s
czestotliwosc_probkowania = 192000;  % Częstotliwość próbkowania w Hz - l. próbek na s
faza_poczotkowa = 0;  % Faza początkowa w radianach
liczba_iteracji = 45000;  % Liczba iteracji
omega = 2 * pi * czestotliwosc / czestotliwosc_probkowania;  % Omega (częstotliwość kątowa)
a = 2 * cos(omega);  % Współczynnik a

y = [sin(faza_poczotkowa + omega), sin(faza_poczotkowa + 2 * omega)];  % Pierwsze dwie próbki
x = zeros(1, liczba_iteracji);  % Tablica na wartości kątów, poczatek w 0
x(2) = omega;  % Druga próbka dla osi X (w radianach)

% Generowanie sygnału rekurencyjnie
for i = 3:liczba_iteracji
   y(i) = a * y(i-1) - y(i-2);  % Rekurencyjna formuła
   x(i) = x(i-1) + omega;  % Zwiększanie kąta, bo omega to częst. kątowa), 
   % wartość, o którą zmienia się argument funkcji sinus przy każdej kolejnej próbce w czasie
end

% Wykres
plot(x, y);
title("Wykres");
xlabel("[rad]");
ylabel("sin(x)");
grid on;


