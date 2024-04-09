# Zmodyfikuj kod programu z poprzedniego zadania w taki sposób, aby wykorzystywał tablice. Sprawdź
# czas działania obu programów (z 1. i 2. zadania). W tym celu wykorzystaj metodę time() z modułu
# time (przykład zostanie wytłumaczony na tablicy).
from array import *
import statistics
import time

ttime = time.time()

L = array('f', [1, 2])

for i in range(2, 48):
    new_el = (L[i-1] + L[i-2]) / (L[i-1] - L[i-2])
    L.append(new_el)
print(L)
print(list(L))

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

print(time.time() - ttime)
