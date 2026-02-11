using JuMP, Ipopt, Plots

# --- Paramètres ---
w = 0.3
x0, y0, θ0 = 0.0, 0.0, π/7
xf, yf, θf = 4.0, 7.0, -π/2
N = 1000

# 3 segments égaux
block_size = ceil(Int, N/3)          # taille approximative d’un segment
segment_indices = [1:block_size, 
                   block_size+1:2*block_size, 
                   2*block_size+1:N]  # le dernier segment ajuste le reste

# --- Modèle JuMP ---
sys = Model(Ipopt.Optimizer)
set_optimizer_attribute(sys, "tol", 1e-6)

@variables(sys, begin
    x[1:N+1]
    y[1:N+1]
    θ[1:N+1]
    -1 ≤ u[1:N] ≤ 1           # contrôles discrets
    -1 ≤ u_block[1:3] ≤ 1     # 3 valeurs constantes par segments
    0 ≤ Δt ≤ 1
end)

@objective(sys, Min, Δt)

# Contraintes de bord
@constraints(sys, begin
    x[1] == x0
    y[1] == y0
    θ[1] == θ0
    x[N+1] == xf
    y[N+1] == yf
    θ[N+1] == θf
end)

# Lier u[j] aux 3 valeurs constantes
for k in 1:3
    for j in segment_indices[k]
        @constraint(sys, u[j] == u_block[k])
    end
end

# Dynamique : implicit midpoint
for j in 1:N
    @NLconstraint(sys, x[j+1] == x[j] + Δt * (w + cos((θ[j] + θ[j+1])/2)))
    @NLconstraint(sys, y[j+1] == y[j] + Δt * sin((θ[j] + θ[j+1])/2))
    @NLconstraint(sys, θ[j+1] == θ[j] + Δt * ((u[j] + u[j])/2))
end

# --- Résolution ---
println("Solving...")
optimize!(sys)

Δt_val = value(Δt)
t = (0:N) * Δt_val
x_sol = value.(x)
y_sol = value.(y)
θ_sol = value.(θ)
u_sol = value.(u_block)

println("u_block optimaux : ", u_sol)
println("Temps final : ", Δt_val * N)

# Préparer vecteur de contrôle pour plot (3 segments égaux)
u_plot_vals = vcat(fill(u_sol[1], length(segment_indices[1])),
                   fill(u_sol[2], length(segment_indices[2])),
                   fill(u_sol[3], length(segment_indices[3])))

# Plots : états
p1 = plot(t, x_sol, xlabel="t", ylabel="x(t)", legend=false)
p2 = plot(t, y_sol, xlabel="t", ylabel="y(t)", legend=false)
p3 = plot(t, θ_sol, xlabel="t", ylabel="θ(t)", legend=false)
p4 = plot(t[1:end-1], u_plot_vals, xlabel="t", ylabel="u(t)", legend=false)
display(plot(p1, p2, p3, p4, layout=(2,2), size=(800,600)))

# Trajectoire
traj = plot(x_sol, y_sol, lw=3, c=:black, xlabel="x", ylabel="y", legend=false, ratio=1)
for i in 1:5:N+1
    scatter!([x_sol[i]], [y_sol[i]], color=:red)
end
display(traj)
