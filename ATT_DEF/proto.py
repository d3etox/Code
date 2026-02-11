import wx
import random

# ----------- Fonctions utilitaires -----------
def tirer_carte():
    return random.randint(1, 13)

# ----------- Classes de jeu -----------
class Player:
    def __init__(self, name, deck):
        self.name = name
        self.hp = 0
        self.shield = 0
        self.charges = []  # cartes face cachée
        self.alive = True
        self.deck = deck    
        self.setup()

    def setup(self):
        cartes = sorted([self.deck.draw() for _ in range(3)])
        self.shield = cartes[0]
        self.hp = cartes[1] + cartes[2]
        self.deck.add_to_board(cartes)

    def total_charges(self):
        return len(self.charges)

    def attack(self, target, log_callback):
        attack_card = self.deck.draw()
        # On révèle les cartes de charge (face cachée)
        revealed_charges = sum(self.charges)
        self.charges.clear()
        damage = max(0, attack_card + revealed_charges - target.shield)
        log_callback(f"{self.name} attaque {target.name} (carte: {attack_card}, charges: {revealed_charges}) → dégâts {damage}")
        if damage > 0:
            target.hp -= damage
            log_callback(f"{target.name} perd {damage} PV → reste {target.hp} PV")
            if target.hp <= 0:
                target.alive = False
                log_callback(f"💀 {target.name} est éliminé !")

    def change_shield(self, target, log_callback):
        new_card = self.deck.draw()
        log_callback(f"{self.name} change le bouclier de {target.name} ({target.shield} → {new_card})")
        target.shield = new_card

    def charge(self, log_callback):
        card = self.deck.draw()
        self.charges.append(card)
        log_callback(f"{self.name} charge (carte cachée) → total {len(self.charges)} charge(s)")

class Deck:
    def __init__(self):
        self.all_cards = [v for v in range(1, 14)] * 4  # 52 cartes
        self.cards = self.all_cards.copy()
        random.shuffle(self.cards)
        self.cards_on_board = []  # cartes actives (ex : PV, bouclier...)

    def draw(self):
        """Tire une carte sans remise, remélange si vide."""
        if not self.cards:
            self.reshuffle()
        card = self.cards.pop()
        return card

    def add_to_board(self, cards):
        """Ajoute des cartes visibles (bouclier, PV, etc.) à la zone de plateau."""
        if isinstance(cards, int):
            cards = [cards]
        self.cards_on_board.extend(cards)

    def remove_from_board(self, cards):
        """Quand une carte quitte le plateau (par ex destruction), on la libère."""
        if isinstance(cards, int):
            cards = [cards]
        for c in cards:
            if c in self.cards_on_board:
                self.cards_on_board.remove(c)

    def reshuffle(self):
        """Remet toutes les cartes non visibles dans la pioche."""
        self.cards = [c for c in self.all_cards if c not in self.cards_on_board]
        random.shuffle(self.cards)
        print("♻️ Nouveau mélange du paquet !")


# ----------- Classe principale du jeu -----------
class Game(wx.Frame):
    def __init__(self, num_players=4):
        super().__init__(None, title="Jeu de cartes - Prototype", size=(800, 600))

        self.deck = Deck()
        self.players = [Player("Charles", self.deck)] + [Player(f"Joueur {i+1}", self.deck) for i in range(1, num_players)]
        self.human_player = self.players[0] 


        # Interface
        panel = wx.Panel(self)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Tableau des joueurs
        self.stats_panel = wx.Panel(panel)
        stats_sizer = wx.GridSizer(rows=num_players, cols=1, hgap=5, vgap=5)
        self.player_labels = []

        for p in self.players:
            lbl = wx.StaticText(self.stats_panel, label=f"{p.name} → PV={p.hp} | Bouclier={p.shield} | Charges={len(p.charges)}")
            self.player_labels.append(lbl)
            stats_sizer.Add(lbl, 0, wx.LEFT, 10)

        self.stats_panel.SetSizer(stats_sizer)
        main_sizer.Add(self.stats_panel, 0, wx.ALL | wx.EXPAND, 10)

        # Zone de logs
        self.log = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY)
        main_sizer.Add(self.log, 1, wx.EXPAND | wx.ALL, 10)

        # Zone info joueur courant
        self.status = wx.StaticText(panel, label="")
        main_sizer.Add(self.status, 0, wx.LEFT | wx.TOP, 10)

        # Boutons d'action
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_attack = wx.Button(panel, label="⚔️ Attaquer")
        self.btn_shield = wx.Button(panel, label="🛡️ Changer un bouclier")
        self.btn_charge = wx.Button(panel, label="⚡ Charger")

        for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
            btn_sizer.Add(btn, 0, wx.ALL, 5)
            btn.Disable()  # désactivé avant le début

        main_sizer.Add(btn_sizer, 0, wx.CENTER)

        panel.SetSizer(main_sizer)

        # Bindings
        self.btn_attack.Bind(wx.EVT_BUTTON, self.on_attack)
        self.btn_shield.Bind(wx.EVT_BUTTON, self.on_shield)
        self.btn_charge.Bind(wx.EVT_BUTTON, self.on_charge)

        # Setup initial
        self.log_message("🃏 Début de la partie !")
        for p in self.players:
            self.log_message(f"{p.name} → PV={p.hp}, Bouclier={p.shield}, Charges={len(p.charges)}")


        self.sort_players()
        self.turn_index = self.players.index(self.human_player)
        self.current_player = self.players[self.turn_index]
        self.start_turn()

        self.Show()

    def update_player_labels(self):
        for i, p in enumerate(self.players):
            status = f"{p.name} → PV={p.hp} | Bouclier={p.shield} | Charges={len(p.charges)}"
            if not p.alive:
                status += " 💀"
            self.player_labels[i].SetLabel(status)

    # ---------- Fonctions de jeu ----------
    def log_message(self, msg):
        self.log.AppendText(msg + "\n")

    def sort_players(self):
        self.players.sort(key=lambda p: (p.hp, p.shield))

    def get_alive_players(self):
        return [p for p in self.players if p.alive]

    def next_turn(self):
        # Vérifie victoire
        alive = self.get_alive_players()
        if len(alive) <= 1:
            winner = alive[0]
            self.log_message(f"\n🏆 {winner.name} remporte la partie avec {winner.hp} PV !")
            for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
                btn.Disable()
            return

        # Avance dans la liste complète des joueurs
        while True:
            self.turn_index = (self.turn_index + 1) % len(self.players)
            self.current_player = self.players[self.turn_index]
            if self.current_player.alive:
                break

        self.start_turn()

    def start_turn(self):
        p = self.current_player
        if not p.alive:
            self.next_turn()
            return
        self.status.SetLabel(f"Tour de {p.name} — PV={p.hp} | Bouclier={p.shield} | Charges={len(p.charges)}")
        self.log_message(f"\n--- Tour de {p.name} ---")

        if p == self.human_player:
            for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
                btn.Enable()
        else:
            wx.CallLater(1000, self.ai_turn, p)

    def ai_turn(self, player):
        if not player.alive:
            self.next_turn()
            return
        action = random.choice(["attack", "shield", "charge"])
        targets = [p for p in self.get_alive_players() if p != player]
        if not targets:
            self.next_turn()
            return
        target = random.choice(targets)

        if action == "attack":
            player.attack(target, self.log_message)
            player.change_shield(target, self.log_message)
        else:
            player.charge(self.log_message)
        self.update_player_labels()
        wx.CallLater(1000, self.next_turn)

    # ---------- Actions joueur humain ----------
    def on_attack(self, event):
        targets = [p for p in self.get_alive_players() if p != self.current_player]
        if not targets:
            return
        choices = [p.name for p in targets]
        dlg = wx.SingleChoiceDialog(self, "Choisir une cible :", "Attaque", choices)
        if dlg.ShowModal() == wx.ID_OK:
            idx = dlg.GetSelection()
            target = targets[idx]
            self.current_player.attack(target, self.log_message)
            for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
                btn.Disable()
            self.update_player_labels()
            wx.CallLater(1000, self.next_turn)
        else:
            # ❌ Si on annule, on ne passe pas le tour : on réactive les boutons
            self.log_message(f"{self.current_player.name} annule l'attaque et réfléchit encore...")
            for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
                btn.Enable()
        dlg.Destroy()


    def on_shield(self, event):
        targets = [p for p in self.get_alive_players()]
        choices = [p.name for p in targets]
        dlg = wx.SingleChoiceDialog(self, "Changer le bouclier de :", "Bouclier", choices)
        if dlg.ShowModal() == wx.ID_OK:
            idx = dlg.GetSelection()
            target = targets[idx]
            self.current_player.change_shield(target, self.log_message)
            for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
                btn.Disable()
            self.update_player_labels()
            wx.CallLater(1000, self.next_turn)
        else:
            # ❌ Si on annule, on reste sur notre tour
            self.log_message(f"{self.current_player.name} renonce à changer de bouclier pour l’instant.")
            for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
                btn.Enable()
        dlg.Destroy()


    def on_charge(self, event):
        self.current_player.charge(self.log_message)
        for btn in [self.btn_attack, self.btn_shield, self.btn_charge]:
            btn.Disable()
        self.update_player_labels()
        wx.CallLater(1000, self.next_turn)

# ----------- Lancement -----------
if __name__ == "__main__":
    app = wx.App()
    Game(num_players=4)
    app.MainLoop()
