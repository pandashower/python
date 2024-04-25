import sys
import time

start = time.time()
moves = []


class Stack:
    def __init__(self, capacity):
        self.capacity = capacity
        self.top = -1
        self.array = [0] * capacity


def createStack(capacity):
    stack = Stack(capacity)
    return stack


def isFull(stack):
    return (stack.top == (stack.capacity - 1))


def isEmpty(stack):
    return (stack.top == -1)


def push(stack, item):
    if (isFull(stack)):
        return
    stack.top += 1
    stack.array[stack.top] = item


def Pop(stack):
    if (isEmpty(stack)):
        return -sys.maxsize
    Top = stack.top
    stack.top -= 1
    return stack.array[Top]


def moveDisksBetweenTwoPoles(src, dest, s, d):
    pole1TopDisk = Pop(src)
    pole2TopDisk = Pop(dest)

    global move_counter
    move_counter += 1

    # kiedy pierwszy wolny
    if (pole1TopDisk == -sys.maxsize):
        push(src, pole2TopDisk)
        moveDisk(d, s)

    # kiedy 2 wolny
    elif (pole2TopDisk == -sys.maxsize):
        push(dest, pole1TopDisk)
        moveDisk(s, d)

    elif (pole1TopDisk > pole2TopDisk):
        push(src, pole1TopDisk)
        push(src, pole2TopDisk)
        moveDisk(d, s)

    else:
        push(dest, pole2TopDisk)
        push(dest, pole1TopDisk)
        moveDisk(s, d)


def moveDisk(fromPeg, toPeg):
    moves.append((fromPeg, toPeg))
    print("Przesuwamy krążek z", fromPeg, "do", toPeg,)


def tohIterative(num_of_disks, src, aux, dest):
    s, d, a = 'A', 'C', 'B'

    # jesli liczba dysków parzysta podmieniamy
    if (num_of_disks % 2 == 0):
        temp = d
        d = a
        a = temp

    for i in range(num_of_disks, 0, -1):
        push(src, i)

    total_num_of_moves = int(pow(2, num_of_disks) - 1)
    for i in range(1, total_num_of_moves + 1):
        if (i % 3 == 1):
            moveDisksBetweenTwoPoles(src, dest, s, d)

        elif (i % 3 == 2):
            moveDisksBetweenTwoPoles(src, aux, s, a)

        elif (i % 3 == 0):
            moveDisksBetweenTwoPoles(aux, dest, a, d)
    return total_num_of_moves


num_of_disks = 25  # robiło się kilka godzin chyba i nic dla 30

src = createStack(num_of_disks)
dest = createStack(num_of_disks)
aux = createStack(num_of_disks)

move_counter = 0

total_theo = tohIterative(num_of_disks, src, aux, dest)

end = time.time()

print("Dla podejścia iteracyjnego gdzie A-src, B-buff, C-dst")
print("Czas wykonania:", end-start)
print("Liczba ruchów:", move_counter, ", Teretyczna optymalna liczba ruchów:", total_theo)
# print(moves)
