import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

ROWS, COLS = 6, 7
PLAYER = 1
AI = 2
DIFFICULTY_LEVELS = {'Easy': 2, 'Medium': 4, 'Hard': 6, 'Impossible': 8}
difficulty = DIFFICULTY_LEVELS['Medium']

board = np.zeros((ROWS, COLS), dtype=int)
game_over = False
turn = PLAYER

fig = plt.figure(figsize=(7,7), facecolor='black')
gs = fig.add_gridspec(2, 1, height_ratios=[6, 1])
ax = fig.add_subplot(gs[0], facecolor='black')
ax_prob = fig.add_subplot(gs[1], facecolor='black')
plt.ion()

# --- Dessin du plateau ---
def draw_board(board, highlight_col=None):
    ax.clear()
    ax.add_patch(Rectangle((-0.5, -0.5), COLS, ROWS, facecolor='blue'))
    for r in range(ROWS):
        for c in range(COLS):
            color = 'black'
            if board[r, c] == PLAYER:
                color = 'yellow'
            elif board[r, c] == AI:
                color = 'red'
            y = ROWS - 1 - r
            circle = Circle((c, y), 0.4, color=color, ec='black', lw=2)
            ax.add_patch(circle)
    if highlight_col is not None:
        rect = Rectangle((highlight_col-0.5, -0.5), 1, ROWS, facecolor='gray', alpha=0.3)
        ax.add_patch(rect)
    ax.set_xlim(-0.5, COLS-0.5)
    ax.set_ylim(-0.5, ROWS-0.5)
    ax.set_xticks(range(COLS))
    ax.set_yticks(range(ROWS))
    ax.set_aspect('equal')
    plt.draw()

# --- Vérification victoire ---
def winning_move(board, piece):
    for c in range(COLS - 3):
        for r in range(ROWS):
            if all(board[r, c+i] == piece for i in range(4)):
                return True
    for c in range(COLS):
        for r in range(ROWS - 3):
            if all(board[r+i, c] == piece for i in range(4)):
                return True
    for c in range(COLS - 3):
        for r in range(ROWS - 3):
            if all(board[r+i, c+i] == piece for i in range(4)):
                return True
        for r in range(3, ROWS):
            if all(board[r-i, c+i] == piece for i in range(4)):
                return True
    return False

def get_valid_locations(board):
    return [c for c in range(COLS) if board[0, c]==0]

def next_open_row(board, col):
    for r in reversed(range(ROWS)):
        if board[r, col]==0:
            return r

def update_probability(score):
    # Normalisation du score avec sigmoïde
    prob_ai = 1 / (1 + np.exp(-score / 200))  # entre 0 et 1

    ax_prob.clear()
    ax_prob.set_xlim(0, 1)
    ax_prob.set_ylim(0, 1)
    ax_prob.axis('off')

    # Barres de probabilité
    ax_prob.barh(0.5, prob_ai, color='red')
    ax_prob.barh(0.5, 1-prob_ai, left=prob_ai, color='yellow')

    # Texte centré
    ax_prob.text(0.02, 0.5, f"{int((1-prob_ai)*100)}% 🟡", va='center', ha='left', color='white', fontsize=12)
    ax_prob.text(0.98, 0.5, f"🔴 {int(prob_ai*100)}%", va='center', ha='right', color='white', fontsize=12)

    plt.draw()


# --- Heuristique améliorée ---
def score_window(window, piece):
    score = 0   
    opp_piece = PLAYER if piece==AI else AI
    if window.count(piece)==4:
        score += 1000
    elif window.count(piece)==3 and window.count(0)==1:
        score += 10
    elif window.count(piece)==2 and window.count(0)==2:
        score += 5
    if window.count(opp_piece)==3 and window.count(0)==1:
        score -= 80  # bloquer l'adversaire
    return score

def score_position(board, piece):
    score = 0
    # Centre
    center_array = list(board[:, COLS//2])
    score += center_array.count(piece)*6
    # Horizontal
    for r in range(ROWS):
        row_array = list(board[r, :])
        for c in range(COLS-3):
            window = row_array[c:c+4]
            score += score_window(window, piece)
    # Vertical
    for c in range(COLS):
        col_array = list(board[:, c])
        for r in range(ROWS-3):
            window = col_array[r:r+4]
            score += score_window(window, piece)
    # Diagonales
    for r in range(ROWS-3):
        for c in range(COLS-3):
            window = [board[r+i, c+i] for i in range(4)]
            score += score_window(window, piece)
        for c in range(COLS-3):
            window = [board[r+3-i, c+i] for i in range(4)]
            score += score_window(window, piece)
    return score

# --- Minimax alpha-bêta ---
def minimax(board, depth, alpha, beta, maximizingPlayer):
    valid_locations = get_valid_locations(board)
    is_terminal = winning_move(board, PLAYER) or winning_move(board, AI) or len(valid_locations)==0
    if depth==0 or is_terminal:
        if winning_move(board, AI):
            return (None, 100000)
        elif winning_move(board, PLAYER):
            return (None, -100000)
        else:
            return (None, score_position(board, AI))
    if maximizingPlayer:
        value = -np.inf
        best_col = np.random.choice(valid_locations)
        for col in valid_locations:
            row = next_open_row(board, col)
            temp_board = board.copy()
            temp_board[row, col] = AI
            new_score = minimax(temp_board, depth-1, alpha, beta, False)[1]
            if new_score > value:
                value = new_score
                best_col = col
            alpha = max(alpha, value)
            if alpha >= beta:
                break
        return best_col, value
    else:
        value = np.inf
        best_col = np.random.choice(valid_locations)
        for col in valid_locations:
            row = next_open_row(board, col)
            temp_board = board.copy()
            temp_board[row, col] = PLAYER
            new_score = minimax(temp_board, depth-1, alpha, beta, True)[1]
            if new_score < value:
                value = new_score
                best_col = col
            beta = min(beta, value)
            if alpha >= beta:
                break
        return best_col, value

# --- Animation chute jeton ---
def animate_drop(piece, col, row):
    for r in range(row+1):
        draw_board(board)
        temp_board = board.copy()
        temp_board[r, col] = piece
        draw_board(temp_board)
        plt.pause(0.05)

# --- Gestion clic ---
def on_click(event):
    global board, turn, game_over
    if game_over or turn != PLAYER or event.xdata is None:
        return
    col = int(event.xdata)
    if board[0, col]!=0:
        return
    row = next_open_row(board, col)
    animate_drop(PLAYER, col, row)
    board[row, col] = PLAYER
    if winning_move(board, PLAYER):
        draw_board(board)
        print("Vous gagnez !")
        game_over = True
        return
    turn = AI
    draw_board(board)
    # IA joue
    col_ai, score_ai = minimax(board, difficulty, -np.inf, np.inf, True)
    update_probability(score_ai)
    if col_ai is not None:
        row_ai = next_open_row(board, col_ai)
        animate_drop(AI, col_ai, row_ai)
        board[row_ai, col_ai] = AI
        if winning_move(board, AI):
            draw_board(board)
            print("IA gagne !")
            game_over = True
            return
    turn = PLAYER
    draw_board(board)

# --- Highlight colonne ---
highlight_col = None
def on_motion(event):
    global highlight_col
    if game_over or event.xdata is None:
        highlight_col = None
        draw_board(board)
        return
    highlight_col = int(event.xdata)
    draw_board(board, highlight_col)

# --- Initialisation ---
draw_board(board)
update_probability(0)
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('motion_notify_event', on_motion)
plt.show(block=True)
