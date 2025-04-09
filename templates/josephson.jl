cd(@__DIR__)
using Pkg; Pkg.activate(".")
using OrdinaryDiffEq,Plots,LaTeXStrings

#######################
# Blatt 2 -- Exercise 3
#######################

# EoM of Josephson junction
# u[1] = ϕ
# p[1] = I/I_c
function josephson!(du,u,p,t)
    du[1]= (p[1]-sin(u[1]))
    return nothing
end

function josephson(u,p,t)
    return (p[1]-sin(u[1]))
end

# define initial position, parameters, and timespan
u0 = [0.0f0]
tspan = (0.0f0, 100.0f0)
p= [1.1f0]

# solve and plot ϕ
prob = ODEProblem(josephson!,u0,tspan,p)
sol = solve(prob,Tsit5())
plot(sol,xaxis=(L"t"),label=L"\phi(t)"*"; "*L"I/I_c="*"$(p[1])",legend=:topleft,lw=2,dpi=300)

# calculate supercurrent and plot on top
supercurrent = map(t->sin(sol(t)[1]),0.0f0:0.01f0:100.0f0)
plot!(collect(0.0f0:0.01f0:100.0f0),supercurrent,xaxis=(L"t"),label=L"\sin(\phi(t))"*"; "*L"I/I_c="*"$(p[1])",legend=:topleft,lw=2,dpi=300)

# also plot ̇ϕ on top 
ϕ = map(t->josephson(sol(t),p,0.0f0),0.0f0:0.01f0:100.0f0)
plot!(collect(0.0f0:0.01f0:100.0f0),ϕ,xaxis=(L"t"),label=L"\dot{\phi}(t)"*"; "*L"I/I_c="*"$(p[1])",legend=:topleft,lw=2,dpi=300)

savefig("./jj_sketch_$(p[1]).png")

# let's also look at the average voltage

# exact expression for the average voltage
# set R = I_c = 1
function V_avg_exact(I)
    return sqrt(I^2-1)
end

# average voltage obtained from numerics 
V_avg_num = sum(ϕ)/length(ϕ)

# exact average voltage
V_avg_ex = V_avg_exact(p[1])

# let's check for a range of values of I/I_c
I_range = collect(1.0f0:0.1f0:100.0f0)

function V_avg_numerical(prob,I)
    p = [I]
    _prob = remake(prob,p=p)
    sol = solve(_prob,Tsit5())
    ϕ = map(t->josephson(sol(t),p,0.0f0),0.0f0:0.01f0:100.0f0)
    V_avg_num = sum(ϕ)/length(ϕ)

    return V_avg_num
end

V_avg_num = map(I -> V_avg_numerical(prob,I), I_range)
V_avg_ex = map(I -> V_avg_exact(I), I_range)

plot(I_range,V_avg_num,xaxis=(L"I/I_c"),yaxis=(L"\langle V/RI_{c} \rangle"),label="numerical",legend=:topleft,lw=2,dpi=300)
plot!(I_range,V_avg_ex,xaxis=(L"I/I_c"),yaxis=(L"\langle V/RI_{c} \rangle"),label="analytical",legend=:topleft,lw=2,dpi=300)
savefig("./jj_IV_curve.png")

# let's zoom in a bit for smaller values of I/I_c
I_range = collect(0.0f0:0.01f0:2.0f0)
I_range_ex = collect(1.0f0:0.01f0:2.0f0)

function V_avg_numerical(prob,I,tstart,dt,tmax)
    _prob = remake(prob,p=[I],tspan=(0.0f0,tmax))
    _sol = solve(_prob,Tsit5(),reltol=1e-9)
    ϕ = map(t->josephson(_sol(t),[I,1.0f0],t),tstart:dt:tmax)
    V_avg_num = sum(ϕ)/length(ϕ)

    return V_avg_num
end

tstart = 100.0f0
dt = 0.01f0
tmax = 1000.0f0

V_avg_num = map(I -> V_avg_numerical(prob,I,tstart,dt,tmax), I_range)
V_avg_ex = map(I -> V_avg_exact(I), I_range_ex)

plot(I_range,V_avg_num,xaxis=(L"I/I_c"),yaxis=(L"\langle V/RI_{c} \rangle"),label="numerical",legend=:topleft,lw=2,dpi=300)
plot!(I_range_ex,V_avg_ex,xaxis=(L"I/I_c"),yaxis=(L"\langle V/RI_{c} \rangle"),label="analytical",legend=:topleft,lw=2,dpi=300)
savefig("./jj_IV_curve_2.png")
