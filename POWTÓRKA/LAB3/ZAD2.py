# Napisz wªasny moduª w j¦zyku Python, który:
# 1.
# b¦dzie zawieraª deklaracj¦ trzech klas: circle, triangle, square,
# 2.
# dla ka»dej z klas b¦dzie posiadaª metody licz¡ce pole powierzchni i obwód gur.
# Zaimportuj napisany moduª. Przetestuj jego dziaªanie.

import mojmodul as mm

kolko = mm.Circle(2)
kolko.ppowierzchni()
kolko.obwod()

kwadrat = mm.Square(5)
kwadrat.ppowierzchni()
kwadrat.obwod()

trojkat = mm.Triangle(3, 4, 5)
trojkat.ppowierzchni()
trojkat.obwod()
