# Napisz dowolny kod, który wywołuje wyjątki:
# IndexError
# ZeroDivisionError
# NameError
# a następnie obsłuż je za pomocą bloków try-except, aby program nie zakończył się z ich powodu

try:
    b = [1, 2]
    print(b[3])
except IndexError as e:
    print(f"IndexError chuju: {e}")

try:
    a = 1/0
except ZeroDivisionError as e1:
    print(f"zero divison miernoto: {e1}")


try:
    print(zmienna_ktorej_nie_ma)
except NameError as e2:
    print(f"WTF nie ma zmiennej?! ale gówno: {e2}")