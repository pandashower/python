# 1.
# wczyta list¦ plików w katalogu 'zadanie1',
# 2.
# utworzy nowe katalogi i przeniesie do nich pliki, których nazwa rozpoczyna si¦ na t¦ sam¡ liter¦
# (np. plik 'eHszo' zostanie przeniesiony do katalogu 'E').
# Podczas realizacji zadania 1 przydatne b¦d¡ takie metody jak: glob.glob(), os.mkdir, os.rename().
# W zadaniu u»yj obsªugi wyj¡tków zamiast sprawdzania czy katalog istnieje
import os
import glob
import shutil

# os.makedirs("./zaddddd/okk")
files = glob.glob(r".\zadanie1\*")

print(files)

for filepath in files:
    first_letter = os.path.basename(filepath)[0]
    newdir = os.path.join(r".\zadanie1", first_letter)

    try:
        os.mkdir(newdir)
    except FileExistsError:
        pass

    newpath = os.path.join(newdir, os.path.basename(filepath))
    shutil.move(filepath, newpath)


