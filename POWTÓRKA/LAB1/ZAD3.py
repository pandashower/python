# Na wykładzie przedstawiono kilka sposobów konstruowaniu pętli ’for’ w języku Python.
# Napisz kod programu który:
# 1.
# implementuje pętle ’for’ na obiekcie iterowanym
# 2.
# implementuje pętle ’for’ w sposób zaczerpnięty z C++
# 3.
# sprawdź czy obie implementacje zwracają takie same wyniki
# Sprawdź czas wykonywania obu pętli.
import time

L = []
for i in range(0, 1000):
    L.append(i)

ptime = time.time()
for items in L:
    L[items] += 1
    print(items)
petime = time.time()

ctime = time.time()
for i in range(0, len(L)):
    L[i] += 1
    print(L[i])
cetime = time.time()

print(f"""Czas for z pythona: {petime - ptime}, 
Czas z c++: {cetime - ctime}""")

