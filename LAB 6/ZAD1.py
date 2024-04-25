import time

start = time.time()


def Hanoi_re(n, sour, dest, buff, moves=[]):
    if n == 1:
        moves.append((sour, dest))
        print(f"Przesuwamy krążek z {sour} do {dest}")
        return 1
    else:
        count = 0
        count += Hanoi_re(n-1, sour, buff, dest, moves)
        moves.append((sour, dest))
        print(f"Przesuwamy krążek z {sour} do {dest}")
        count += 1
        count += Hanoi_re(n-1, buff, dest, sour, moves)
        return count

# Przykładowe użycie
number_of_disks = 25  # 8 do testu
moves = []
total_steps = Hanoi_re(number_of_disks, 'A', 'C', 'B', moves)
total_theoretical = int(pow(2, number_of_disks) - 1)

end = time.time()

print("Dla podejścia rekurencyjnego gdzie A-src, B-buff, C-dst")
print("Czas wykonania:", end-start)
print(f"Liczba ruchów: {total_steps}", ", Teretyczna optymalna liczba ruchów:", total_theoretical)
# print("Lista ruchów:", moves)
