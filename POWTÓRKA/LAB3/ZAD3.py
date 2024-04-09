# Caªkowanie metod¡ Monte Carlo polega na:
# 1.
# okre±leniu prostok¡ta w którym znajduje si¦ caªkowany przedziaª,
# 2.
# losowaniu wspóªrz¦dnych wewn¡trz prostok¡ta,
# 3.
# sprawdzeniu czy punkt z wylosowanych wspóªrz¦dnych speªnia równanie,
# 4.
# obliczenia stosunku ilo±ci punktów speªniaj¡cych równanie do wszystkich wylosowanych punk-
# tów,
# 5.
# obliczenie warto±ci caªki na podstawie stosunku ilo±ci traonych punktów do pola pola prosto-
# k¡ta.
# Napisz wªasn¡ implementacj¦ caªkowania wedªug metody Monte Carlo. Przetestuj dziaªanie
# programu obliczaj¡c pole koªa o promieniu R oraz całka od 0 do 2 z sin(x)
# Jak zmienia si¦ dokªadno±¢ oszacowania wraz ze wzrostem ilo±ci wylosowanych punktów?
# Przydatne b¦d¡ metody z takich moduªów jak random oraz math
import random
import math

# sin(x)
c = 0  # licznik punktów pod funkcją
b = 2  # ograniczenie x
g = 1  # górna granica y dla sin(x)
N = 10000  # wielkość próby

for i in range(N):
    x = random.uniform(0, b)
    y = random.uniform(0, g)
    if y < math.sin(x):
        c += 1
print("Pole całki oznaczonej sin(x) od 0 do 2: ", (b*g*c)/N)


# kolo
c = 0
R = int(input("Podaj promień koła: "))
b = R * 2
g = R * 2
N = 10000

for i in range(N):
    x = random.uniform(0, b)
    y = random.uniform(0, g)
    if R < x*x + y*y:
        c += 1
print("Pole Twojego koła: ", (b*g*c)/N)



# def monte_carlo_circle(radius, num_points):
#     inside_circle = 0
#     total_points = num_points
#
#     for _ in range(num_points):
#         x = random.uniform(-radius, radius)
#         y = random.uniform(-radius, radius)
#         if x**2 + y**2 <= radius**2:
#             inside_circle += 1
#
#     # Pole koła = stosunek punktów w kole do wszystkich punktów * pole prostokąta
#     circle_area = (inside_circle / total_points) * (4 * radius**2)
#     return circle_area
#
# def monte_carlo_integral(func, lower_limit, upper_limit, num_points):
#     inside_function = 0
#     total_points = num_points
#
#     for _ in range(num_points):
#         x = random.uniform(lower_limit, upper_limit)
#         y = random.uniform(0, 1)  # Ograniczenie y do [0, 1], ponieważ sin(x) mieści się w zakresie [0, 1]
#         if y <= func(x):
#             inside_function += 1
#
#     # Całka oznaczona = stosunek punktów pod funkcją do wszystkich punktów * pole prostokąta
#     integral_value = (inside_function / total_points) * (upper_limit - lower_limit)
#     return integral_value
#
# # Definicja funkcji sin(x)
# def sine_function(x):
#     return math.sin(x)
#
# # Ustawienia
# radius = 1
# lower_limit = 0
# upper_limit = 2
# num_points = 10000  # Zwiększ tę wartość, aby uzyskać dokładniejsze oszacowanie
#
# # Obliczenia
# circle_area = monte_carlo_circle(radius, num_points)
# integral_value = monte_carlo_integral(sine_function, lower_limit, upper_limit, num_points)
#
# # Wyniki
# print("Pole koła o promieniu", radius, "to:", circle_area)
# print("Całka oznaczona sin(x) w przedziale od", lower_limit, "do", upper_limit, "to:", integral_value)
