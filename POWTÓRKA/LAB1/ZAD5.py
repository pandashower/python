# Napisz kod gry w kółko i krzyżyk dla dwóch graczy na planszy 3x3
# 1.
# Napisz funkcję, która wyświetla planszę w postaci:
# -------------
# | | | |
# -------------
# | | | |
# -------------
# | | | |
# -------------
# 2.
# Napisz funkcję, która wczyta ruchy graczy i je wyświetli je na ekranie komputera.
# 3.
# Napisz funkcję, która sprawdzi po każdym ruchu czy rozgrywka się zakończyła i wyświetli wynik.
# Zaimplementuj algorytm gry w kółko i krzyżyk tak, aby można było zagrać z komputerem

def plansza(board):
    A, B, C, D, E, F, G, H, J = (board[0][0], board[1][0], board[2][0], board[0][1],
                                 board[1][1], board[2][1], board[0][2], board[1][2], board[2][2])
    print(f"""
# -------------
# | {A} | {B} | {C} |
# -------------
# | {D} | {E} | {F} |
# -------------
# | {G} | {H} | {J} |
# ------------- """)


def check_winner(board, player):
    #poziomo
    if (board[0][0] == player and board[1][0] == player and board[2][0] == player
       or board[0][1] == player and board[1][1] == player and board[2][1] == player
       or board[0][2] == player and board[1][2] == player and board[2][2] == player):
        return True  #pionowo
    elif (board[0][0] == player and board[0][1] == player and board[0][2] == player
          or board[1][0] == player and board[1][1] == player and board[1][2] == player
          or board[2][0] == player and board[2][1] == player and board[2][2] == player):
        return True  #ukos
    elif (board[0][0] == player and board[1][1] == player and board[2][2] == player
          or board[0][2] == player and board[1][1] == player and board[2][0] == player):
        return True

def graj():
    board = [[" ", " ", " "],
             [" ", " ", " "],
             [" ", " ", " "]]
    plansza(board)

    current_player = "X"

    for i in range(0, 9):
        while True:
            x = int(input(f"Graczu {current_player} podaj kolumne 0-2: "))
            y = int(input(f"Graczu {current_player} podaj wiersz 0-2: "))
            if 0 <= x < 3 and 0 <= y < 3:
                break
            else:
                print("Niepoprawny ruch!")

        board[x][y] = current_player
        plansza(board)

        if check_winner(board, current_player):
            print(f"Gracz {current_player} wygrywa!")
            return

        if current_player == "X":
            current_player = "O"
        elif current_player == "O":
            current_player = "X"

    print("Remis!")


graj()
