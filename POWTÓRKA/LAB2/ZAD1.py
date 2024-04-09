# Napisz program, który:
# 1.
# znajdzie wszystkie liczby podzielne przez 7, które s¡ jednocze±nie niepodzielne przez 5 z prze-
# dziaªu <500,3000>
# 2.
# zapisze liczby w postaci ci¡gu znakowego (string) bez przerw,
# 3.
# policzy wyst¡pienia ci¡gu '21' i zamieni je na 'XX'

ciag_znakowy = ""

for i in range(500, 3001):
    if i % 7 == 0 and i % 5 != 0:
        ciag_znakowy += str(i)

wystapienia = ciag_znakowy.count('21')
ciag_znakowy = ciag_znakowy.replace("21", "XX")

print(ciag_znakowy)
print(f"wystąpienia 21: {wystapienia}")

