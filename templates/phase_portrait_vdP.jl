cd(@__DIR__)
using Pkg; Pkg.activate(".")
using OrdinaryDiffEq,Plots,LaTeXStrings

#######################
# Blatt 1 -- Exercise 5
#######################

# equations of motion of the van der Pol oscillator
# u1 = x
# u2 = v
# p = eta
function vdP!(du,u,p,t)
    du[1]= u[2]
    du[2] = -u[1]-p*u[2]*(u[1]^2-1)
    return nothing
end

# set eta parameter 
p = -5.0f0

# set initial condition
x0 = 0.9f0
v0 = 0.0f0
u0 = [x0,v0]

# set grid for iteration of initial conditions
x0min = -3.0f0
x0max = 3.0f0
dx = 1.0f0
x_range = collect(x0min:dx:x0max)

v0min = -4.0f0
v0max = 4.0f0
dv = 1.0f0
v_range = collect(v0min:dv:v0max)

# set timespan
T = 30.0f0
tspan = (0.0f0,T)

# solve equations of motion
prob = ODEProblem(vdP!,u0,tspan,p)
sol = solve(prob,Tsit5())

# plot
plt = plot(sol, idxs=(1, 2), title=L"\eta ="*"$p", xaxis=L"x", yaxis=L"\dot{x}", leg=false, xlims=(x0min-2.0f0,x0max+2.0f0),ylims=(v0min-2.0f0,v0max+2.0f0),dpi=300)

# iterate over different initial conditions
for v0 in v_range
    for x0 in x_range
        # remake ODEProblem with updated initial condition and parameter value
        _prob = remake(prob, u0=[x0,v0], p=p)

        # solve and plot
        sol = solve(_prob, Tsit5())
        plot!(plt, sol, idxs=(1, 2), xlims=(x0min-2.0f0,x0max+2.0f0),ylims=(v0min-2.0f0,v0max+2.0f0))
    end
end

# show final plot and save it in local directory
plot(plt)
savefig("./vdP_phase_plot_eta=$p.png")