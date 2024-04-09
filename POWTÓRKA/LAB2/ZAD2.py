# wczyta plik zadanie2.csv,
# usunie z niego linie w których warto±¢ 'val' jest pustym ci¡giem znakowym,
# posortuje linie w pliku po 'id',
# 4.
# naprawi numeracj¦ 'id' w taki sposób, aby nie byªo duplikatów (w przypadku znalezienia dupli-
# katu nale»y nada¢ numer o jeden wi¦kszy od poprzedniego),
# 5.
# zamie« wszystkie du»e litery na maªe,
# 6.
# usunie wszystkie wyrazy w których dwuliterowy preks skªada si¦ ze znaków, które s¡ obok siebie
# na tablicy ASCII (np. g-h), wy±wietl usuni¦te wyrazy i ich 'id'.
import csv

with open("zadanie2.csv", "r", newline="") as csv_file:
    csv_reader = csv.reader(csv_file)
    dane = []
    for linie in csv_reader:
        if linie[1]:
            dane.append({'id': linie[0], 'val': linie[1]})

    with open("new.csv", "w", newline="") as new_file:
        csv_writer = csv.writer(new_file)
        for row in dane:
            csv_writer.writerow([row['id'], row['val']])


with open("new.csv", "r", newline="") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    posortowane = sorted(csv_reader, key=lambda x: int(x[0]))

with open("newnew.csv", "w", newline="") as file:
    csv_writer = csv.writer(file)
    csv_writer.writerow(header)
    csv_writer.writerows(posortowane)


# with open("newnew.csv", "r", newline="") as file:
#     csv_reader = csv.reader(file)
#     header = next(csv_reader)
#
#     with open("newnew_NEW.csv", "w", newline="") as new_file:
#         csv_writer = csv.writer(new_file)
#         csv_writer.writerow(header)
#         nowe_dane = []
#         counter = 0
#
#         for linie in csv_reader:
#             counter += 1
#             nowe_dane.append(linie)
#             csv_writer.writerow([id: str(linie), val])

# po prostu po koleji + małe litery
with open("newnew.csv", "r", newline="") as csv_file:
    csv_reader = csv.reader(csv_file)
    header = next(csv_reader)
    dane = []
    for index, linie in enumerate(csv_reader, start=1):
        if linie[1]:
            dane.append({'id': str(index), 'val': linie[1]})
    with open("newnew_NEW.csv", "w", newline="") as new_file:
        csv_writer = csv.writer(new_file)
        csv_writer.writerow(header)
        for row in dane:
            csv_writer.writerow([row['id'], row['val'].lower()])
