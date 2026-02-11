using OptimalControl
using NLPModelsIpopt
using Plots

# Parameters
w = 0.8
x0 = 0 
y0 = 0 
θ0 = π/7
xf = 4
yf = 7
θf =-π/2 
N = 100

# Définition du problème de commande optimale
ocp = @def begin
    tf ∈ R, variable 
    t ∈ [0, tf], time                
    X = (x,y,θ) ∈ R³, state                  
    u ∈ R, control                    
    X(0) == [x0,y0,θ0]
    X(tf) == [xf,yf,θf]
    -1 ≤ u(t) ≤ 1                                  
    Ẋ(t) == [w+cos(θ(t)),sin(θ(t)), u(t)]         
    tf → min                
end

# Résolution du problème avec Ipopt
sol = solve(ocp, disc_method=:gauss_legendre_2)

# Tracé des résultats
plot(sol)