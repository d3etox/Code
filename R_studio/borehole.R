
# ============================================================
# TP – Évaluation du cours de Bertrand Iooss
# Polytech Nice Sophia – 2026
# Modèle de forage BOREHOLE
# Auteur : Charles Fabre
# ============================================================
pdf("R_studio/figures.pdf")

rm(list = ls())
set.seed(123)


# ============================================================
# 0) MODÈLE BOREHOLE (DONNÉ PAR LE PROF — À NE PAS MODIFIER)
# ============================================================

borehole <- function(x)
{
xx <- matrix(x, ncol=13)

rw <- xx[,1]
riw <- xx[,2]
r <- xx[,3]
Tu <- xx[,4]
Hu <- xx[,5]
Tum <- xx[,6]
Hum <- xx[,7]
Tlm <- xx[,8]
Hlm <- xx[,9]
Tl <- xx[,10]
Hl <- xx[,11]
L <- xx[,12]
Kw <- xx[,13]

frac1 <- 2 * pi * Tu * (Hu - Hl)
frac2a <- 2 * L * Tu / (log(r / rw) * rw^2 * Kw)
frac2b <- Tu / Tl
frac2 <- log(r / rw) * (1 + frac2a + frac2b)

y <- frac1 / frac2
return(y)
}


# ============================================================
# 1) ÉCHANTILLONNAGE MONTE CARLO
# ============================================================

EchantBorehole <- function(N){

X <- matrix(NA, N, 13)

X[,1] <- rnorm(N, 0.1, 0.015)
X[,2] <- rnorm(N, 0.05, 0.01)
X[,3] <- rlnorm(N, 7.71, 1.0056)
X[,4] <- runif(N, 63100, 116000)
X[,5] <- runif(N, 1000, 1100)
X[,6] <- runif(N, 6310, 11600)
X[,7] <- runif(N, 900, 1000)
X[,8] <- runif(N, 631, 1160)
X[,9] <- runif(N, 800, 900)
X[,10] <- runif(N, 63.1, 116)
X[,11] <- runif(N, 700, 800)
X[,12] <- runif(N, 1120, 1680)
X[,13] <- runif(N, 3000, 12000)

colnames(X) <- c("rw","riw","r","Tu","Hu","Tum","Hum",
"Tlm","Hlm","Tl","Hl","L","Kw")
return(X)
}


# ============================================================
# 2) QUESTION 1 — ÉVALUATIONS MIN / MOY / MAX
# ============================================================


cat("Début : Valeurs minimales (bornes déterministes)\n")

x_min <- c(0.05, 0.02, 100, 63100, 1000, 6310, 900,
631, 800, 63.1, 700, 1120, 3000)

y_min <- borehole(x_min)
cat("y_min =", y_min, "\n")



cat("Début : Valeurs maximales\n")

x_max <- c(0.15, 0.08, 50000, 116000, 1100, 11600, 1000,
1160, 900, 116, 800, 1680, 12000)

y_max <- borehole(x_max)
cat("y_max =", y_max, "\n")



cat("Début : Valeurs moyennes des lois\n")

mean_r <- exp(7.71 + 1.0056^2 / 2) # moyenne d'une lognormale

x_mean <- c(0.10, 0.05, mean_r,
mean(c(63100,116000)),
mean(c(1000,1100)),
mean(c(6310,11600)),
mean(c(900,1000)),
mean(c(631,1160)),
mean(c(800,900)),
mean(c(63.1,116)),
mean(c(700,800)),
mean(c(1120,1680)),
mean(c(3000,12000)))

y_mean <- borehole(x_mean)
cat("y_mean =", y_mean, "\n")


# ============================================================
# 3) QUESTION 2 — PROPAGATION D’INCERTITUDES
# ============================================================



N <- 1000
cat("Début : Monte Carlo (N =", N, ")\n")
X <- EchantBorehole(N)
Y <- borehole(X)

mean_Y <- mean(Y)
var_Y <- var(Y)
cat("Monte Carlo terminé. Moyenne:", mean_Y, ", Variance:", var_Y, "\n")

hist(Y, breaks = 30, probability = TRUE,
main = "Histogramme du débit (Monte Carlo)",
xlab = "Débit (m3/an)")


# --- 2.b Comparaison E[f(X)] vs f(E[X])
# (à commenter dans le rapport : non-linéarité)


B <- 500
cat("Début : Calcul IC du quantile 95% (", B, " répétitions)\n")
q95 <- quantile(Y, 0.95)
q95_rep <- numeric(B)

for (b in 1:B) {
	Yb <- borehole(EchantBorehole(N))
	q95_rep[b] <- quantile(Yb, 0.95)
	if (b %% 50 == 0) cat("  Répétition :", b, "/", B, "\n")
}
cat("IC du quantile 95% terminé.\n")

IC_q95 <- quantile(q95_rep, c(0.025, 0.975))



cat("Début : Estimation P(Y > 250) par Monte Carlo\n")

Ns <- seq(2e6, 2e7, by = 2e6)
p_hat <- numeric(length(Ns))

for (i in seq_along(Ns)) {
	cat("  Taille échantillon :", Ns[i], "... ")
	Ytmp <- borehole(EchantBorehole(Ns[i]))
	p_hat[i] <- mean(Ytmp > 250)
	cat("p_hat =", p_hat[i], "\n")
}
cat("Estimation P(Y > 250) terminée.\n")

plot(Ns, p_hat, type = "b",
xlab = "Taille de l'échantillon",
ylab = "P(Y > 250)",
main = "Convergence Monte Carlo – Événement rare")


# ============================================================
# 4) QUESTION 3 — ANALYSE DE SENSIBILITÉ
# ============================================================


cat("Début : Scatterplots\n")

N <- 500
X <- EchantBorehole(N)
Y <- borehole(X)

par(mfrow = c(4,4))
for (i in 1:13) {
	plot(X[,i], Y,
			 xlab = colnames(X)[i],
			 ylab = "Débit")
	if (i %% 4 == 0) cat("  Scatterplot variable", i, "/ 13\n")
}
cat("Scatterplots terminés.\n")



cat("Début : Calcul indices SRC²\n")

Xn <- scale(X)
Yn <- scale(Y)

lm_mod <- lm(Yn ~ Xn)
SRC2 <- summary(lm_mod)$coefficients[-1,1]^2
SRC2 <- SRC2 / sum(SRC2)
cat("Indices SRC² :\n")
print(SRC2)


# ============================================================
# 5) QUESTION 3.b — INDICES DE SOBOL
# ============================================================


cat("Chargement du package sensitivity\n")
library(sensitivity)

N <- 2000

cat("Début : Calcul indices de Sobol (", N, " échantillons)\n")
sob <- sobolSalt(model = borehole,
				 X1 = EchantBorehole(N),
				 X2 = EchantBorehole(N),
				 nboot = 200)
cat("Indices de Sobol calculés.\n")

print(sob)
cat("Affichage du plot Sobol\n")
plot(sob, choice = 1)


# ============================================================
# 6) QUESTION 4 — MÉTAMODÈLE DE KRIGEAGE
# ============================================================


cat("Chargement du package DiceKriging\n")
library(DiceKriging)

X <- EchantBorehole(500)
Y <- borehole(X)

cat("Début : Apprentissage du krigeage\n")
krig <- km(~., design = X, response = Y)
cat("Krigeage terminé.\n")


# --- Validation (Q²)

Y_pred <- predict(krig, X, type = "UK")$mean
Q2 <- 1 - mean((Y - Y_pred)^2) / var(Y)
cat("Q2 (validation krigeage) =", Q2, "\n")


# --- Sobol sur le métamodèle


cat("Début : Calcul indices de Sobol sur le métamodèle\n")
sob_krig <- sobolSalt(
	model = function(X) predict(krig, X, type = "UK")$mean,
	X1 = EchantBorehole(N),
	X2 = EchantBorehole(N),
	nboot = 200
)
cat("Indices de Sobol (krigeage) calculés.\n")

cat("Affichage du plot Sobol (krigeage)\n")
plot(sob_krig, choice = 1)



cat("\nFin du script R.\n")
cat("Toutes les étapes principales ont été exécutées.\n")

dev.off()