# Napisz program który:
# 1.
# wczyta z klawiatury dowolny tekst
# 2.
# sprawdzi czy wpisany tekst to jeden wyraz
# 3.
# sprowadzi wyraz do zapisu w postaci maªych liter
# 4.
# sprawdzi czy wpisany wyraz to sªowo (zgodne z SJP)
# 5.
# wy±wietli rezultat i czas przetwarzania
# Jako sªownika u»yj zaª¡czonego pliku SJP.txt. Czy mo»na w jaki± ªatwy sposób skróci¢ czas przetwa-
# rzania danych w programie?
import time

dowolny = input("Podaj dowolny tekst: ")

if dowolny.__contains__(" "):
    print("to nie jest jeden wyraz")
else:
    dowolny = dowolny.lower()

file = "./SJP.txt"
with open(file, "r", encoding="utf-8") as SJP:
    slownik = {slowa.strip() for slowa in SJP}  #or [] for a list (now its a set)

start_time = time.time()
if dowolny in slownik:
    print("to slowo")
else:
    print("to gowno")
end_time = time.time()

print(f"Czas wyszukiwania: {end_time-start_time}")
