#include <iostream>
#include <cmath>  // do funkcji trygonometrycznych
#include <vector>

int main() {
    int czestotliwosc = 20;  // Częstotliwość sygnału w Hz
    int czestotliwosc_probkowania = 192000;  // Częstotliwość próbkowania w Hz
    int liczba_iteracji = 45000;  // Liczba iteracji

    // Stała PI (zamiast pi z MATLAB-a)
    const double PI = 3.14159265358979323846;

    // Tu dać double żeby coś z tego było - Zmienna omega w liczbach całkowitych (przybliżona wartość)
    double omega = 2 * PI * czestotliwosc / czestotliwosc_probkowania;  // Omega (częstotliwość kątowa)
    double a = 2 * cos(omega);  // Współczynnik a, wynik w double

    // Inicjalizacja wektora na próbki sygnału y
    std::vector<double> y(liczba_iteracji);
    y[0] = sin(omega);  // Pierwsza próbka
    y[1] = sin(2 * omega);  // Druga próbka

    // Inicjalizacja kąta dla X
    std::vector<double> x(liczba_iteracji);
    x[1] = omega;  // Zapisanie drugiej próbki

    // Generowanie sygnału rekurencyjnie
    for (int i = 2; i < liczba_iteracji; i++) {
        // Obliczanie wartości próbki sygnału y
        y[i] = a * y[i - 1] - y[i - 2];

        // Zwiększanie kąta
        x[i] = x[i - 1] + omega;
    }

    // Wyświetlenie pierwszych kilku wartości y i x, aby zobaczyć wyniki
    std::cout << "Pierwsze kilka wartosci y (sinus):" << std::endl;
    for (int i = 0; i < 10; i++) {
        std::cout << "y[" << i << "] = " << y[i] << std::endl;
    }

    std::cout << "Pierwsze kilka wartosci x (radiany):" << std::endl;
    for (int i = 0; i < 10; i++) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    // Sprawdzanie wartości a dla 20 Hz i 40 Hz (równości z zadania 2.22 i 2.23)
    int omega_20Hz = 2 * PI * 20 / czestotliwosc_probkowania;
    int omega_40Hz = 2 * PI * 40 / czestotliwosc_probkowania;

    //tu normalnie double
    int a_20Hz = static_cast<int>(2 * cos(omega_20Hz) * pow(2, 14));  // Rzutowanie do liczby całkowitej
    int a_40Hz = static_cast<int>(2 * cos(omega_40Hz) * pow(2, 14));  // Rzutowanie do liczby całkowitej

    std::cout << "Wartosc a dla 20 Hz: " << a_20Hz << std::endl;
    std::cout << "Wartosc a dla 40 Hz: " << a_40Hz << std::endl;

    if (a_20Hz == a_40Hz) {
        std::cout << "Wartosci a sa rowne dla 20 Hz i 40 Hz, co jest niepozodanym efektem." << std::endl;
    } else {
        std::cout << "Wartości a są różne dla 20 Hz i 40 Hz." << std::endl;
    }

    return 0;
}
