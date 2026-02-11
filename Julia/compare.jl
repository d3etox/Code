using JuMP, Ipopt, Plots

# --- Paramètres ---
w = 0.3
x0, y0, θ0 = 0.0, 0.0, π/7
xf, yf, θf = 4.0, 7.0, -π/2
N = 1000

# --- Fonction pour résoudre le problème ---
function solve_boat(N; n_segments=0)
    sys = Model(Ipopt.Optimizer)
    set_optimizer_attribute(sys, "tol", 1e-6)

    @variables(sys, begin
        x[1:N+1]
        y[1:N+1]
        θ[1:N+1]
        -1 ≤ u[1:N] ≤ 1
        0 ≤ Δt ≤ 1
    end)

    # Contraintes de bord
    @constraints(sys, begin
        x[1] == x0
        y[1] == y0
        θ[1] == θ0
        x[N+1] == xf
        y[N+1] == yf
        θ[N+1] == θf
    end)

    # Segments constants
    if n_segments > 0
        block_size = ceil(Int, N/n_segments)
        segment_indices = [1:block_size, block_size+1:2*block_size, 2*block_size+1:N]
        @variables(sys, begin
            u_block[1:n_segments]
        end)
        for k in 1:n_segments
            set_lower_bound(u_block[k], -1)
            set_upper_bound(u_block[k], 1)
            for j in segment_indices[k]
                @constraint(sys, u[j] == u_block[k])
            end
        end
    end

    # Dynamique : implicit midpoint
    for j in 1:N
        @NLconstraint(sys, x[j+1] == x[j] + Δt * (w + cos((θ[j] + θ[j+1])/2)))
        @NLconstraint(sys, y[j+1] == y[j] + Δt * sin((θ[j] + θ[j+1])/2))
        @NLconstraint(sys, θ[j+1] == θ[j] + Δt * ((u[j] + u[j])/2))
    end

    # Objectif
    @objective(sys, Min, Δt)

    optimize!(sys)

    Δt_val = value(Δt)
    t = (0:N) * Δt_val
    x_sol = value.(x)
    y_sol = value.(y)
    θ_sol = value.(θ)
    if n_segments > 0
        u_sol = value.(u_block)
        u_plot_vals = vcat(fill(u_sol[1], length(segment_indices[1])),
                           fill(u_sol[2], length(segment_indices[2])),
                           fill(u_sol[3], length(segment_indices[3])))
    else
        u_plot_vals = value.(u)
    end

    return t, x_sol, y_sol, θ_sol, u_plot_vals, Δt_val*N
end

# --- Résoudre les deux cas ---
t_free, x_free, y_free, θ_free, u_free, tf_free = solve_boat(N, n_segments=0)
t_seg, x_seg, y_seg, θ_seg, u_seg, tf_seg = solve_boat(N, n_segments=3)

println("Temps final sans segment : ", tf_free)
println("Temps final 3 segments : ", tf_seg)

# --- Trajectoire comparée ---
traj = plot(x_free, y_free, lw=2, c=:blue, label="u libre")
plot!(x_seg, y_seg, lw=2, c=:red, label="3 segments")
scatter!([x_free[1], x_seg[1]], [y_free[1], y_seg[1]], color=:black, marker=:star, label="départ")
scatter!([x_free[end], x_seg[end]], [y_free[end], y_seg[end]], color=:green, marker=:star, label="arrivée")
xlabel!("x")
ylabel!("y")
title!("Comparaison des trajectoires")
display(traj)

# --- Contrôles comparés ---
ctrl = plot(t_free[1:end-1], u_free, lw=2, c=:blue, label="u libre")
plot!(t_seg[1:end-1], u_seg, lw=2, c=:red, label="3 segments")
xlabel!("t")
ylabel!("u(t)")
title!("Comparaison des contrôles")
display(ctrl)
