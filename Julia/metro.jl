using JuMP, Plots, Ipopt

# --- Parameters ---
q0 = 0.0
v0 = 15
N = 1000
qf = 100
vf = 0.0

# JuMP model
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)
set_optimizer_attribute(model, "tol", 1e-6)
set_optimizer_attribute(model, "max_iter", 1000)

# Variables
@variables(model, begin
    q[1:N+1]
    v[1:N+1]
    -1 <= u[1:N] <= 1
    0 <= Δt <= 1
end)

# Initial/final conditions
@constraints(model, begin
    q[1] == q0
    v[1] == v0
    q[N+1] == qf
    v[N+1] == vf
end)

# Dynamics: Implicit Midpoint
for j in 1:N
    @NLconstraint(model, q[j+1] == q[j] + Δt * (v[j] + v[j+1])/2)
    @NLconstraint(model, v[j+1] == v[j] + Δt * u[j])
end

# Objective: minimize final time
@objective(model, Min, Δt)

# Solve
println("Solving...")
optimize!(model)

q_sol = value.(q)
v_sol = value.(v)
u_sol = value.(u)
Δt_sol = value(Δt)
t = (0:N) * Δt_sol
println("tf = ", N*Δt_sol)

# Plot states
p1 = plot(t, q_sol, xlabel="t", ylabel="q(t)", legend=false)
p2 = plot(t, v_sol, xlabel="t", ylabel="v(t)", legend=false)
p3 = plot(t[1:end-1], u_sol, xlabel="t", ylabel="u(t)", legend=false)
display(plot(p1, p2, p3; layout=(3,1)))
