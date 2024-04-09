# Napisz kod wykonujący następujące operacje:
# 1.
# Wygeneruj listę L zawierającą 48 elementów.
# 2.
# Każdy element listy tworzony jest jako suma dwóch poprzednich wartości podzielona przez ich
# różnicę, przy czym dwa pierwsze elementy to 1 i 2.
# 3.
# Policz średnią oraz modę wartości z listy.
# 4.
# Wypisz na ekran wartości, które pojawiły się w liście więcej niż raz lub poinformuj o ich braku.
import statistics
import time

ltime = time.time()

L = [1, 2]

for i in range(2, 48):
    new_el = (L[i-1] + L[i-2]) / (L[i-1] - L[i-2])
    L.append(new_el)
print(L)

srednia = sum(L) / len(L)
print("Średnia:", srednia)

try:
    print(statistics.mode(L))
except statistics.StatisticsError:
    print("Brak mody")

duplikaty = []
for i in L:
    if L.count(i) > 1:
        duplikaty.append(i)

if duplikaty:
    print(f"Duplikaty: {duplikaty}")
else:
    print("brak duplikatów")

print(time.time() - ltime)
