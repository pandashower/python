# Napisz wªasny moduª w j¦zyku Python, który:
# 1.
# b¦dzie zawieraª deklaracj¦ trzech klas: circle, triangle, square,
# 2.
# dla ka»dej z klas b¦dzie posiadaª metody licz¡ce pole powierzchni i obwód gur.
# Zaimportuj napisany moduª. Przetestuj jego dziaªanie.
import math


class Circle:
    def __init__(self, r):
        self.r = r

    def ppowierzchni(self):
        print("pole:", math.pi*self.r*self.r)

    def obwod(self):
        print("obwod:", 2*math.pi*self.r)


class Square:
    def __init__(self, a):
        self.a = a

    def ppowierzchni(self):
        print("pole:", self.a*self.a)

    def obwod(self):
        print("obwod", 4*self.a)


class Triangle:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def ppowierzchni(self):
        maks = max(self.a, self.b, self.c)
        if self.a + self.b + self.c - maks <= maks:
            print("to nie trojkat")
        else:
            p = (self.a + self.b + self.c)/2
            print("pole:", math.sqrt(p*(p-self.a)*(p-self.b)*(p-self.c)))

    def obwod(self):
        maks = max(self.a, self.b, self.c)
        if self.a + self.b + self.c - maks <= maks:
            print("to nie trojkat")
        else:
            p = (self.a + self.b + self.c)
            print("obwod:", p)
