using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LaTeXStrings

#=
---------------------------------------------------------------------
Parameters
---------------------------------------------------------------------
name            | description                           | default   |
---------------------------------------------------------------------
b               | Birth rate of susceptible             | 0.1, but always multiplied by pop size (1000)
δ               | Death rate                            | 0.1
β               | Transmission rate                     | 0.005
ν               | Recovery rate                         | 0.2
μ               | Mutation rate                         | 0.05
k               | Strain index                          | N/A
a               | Detection scaling from ref. strain    | 1.0
dₘₐₓ            | Maximum detection rate                | 0 
nₖ              | Number of mutation classes            | 3
dₖ = dₘₐₓ e⁻ᵏᵃ  | Scaled detection rate                 | N/A

-------------------------------------------------------------------
Dynamics 
-------------------------------------------------------------------
Ṡ = b - S(δ + β ∑ₖIₖ)
İ₀ = βSI₀ - (ν + μ + δ + dₖ)I₀
İₖ = βSIₖ - (ν + μ + δ + dₖ)Iₖ + μIₖ₋₁  {∀k | 0 < k < n}
İₙ = βSIₙ - (ν + δ + dₖ)Iₙ + μIₙ₋₁ 
Ṙ = ∑ₖ (ν + dₖ)Iₖ₋₁  - δ*R

=# 


#=
    Implement dX/dt 
=#

detectability(dₘₐₓ, a, k::Number) = dₘₐₓ * exp(-a*k)
detectability(dₘₐₓ, a, k::Vector) = dₘₐₓ * exp.(-a.*k)

function SIR(u, θ, t)
    b, δ, β, ν, μ, a, dₘₐₓ, nₖ = θ
    S, I₀, Iₖ, Iₙ, R = u[begin], u[2], u[3:end-2], u[end-1], u[end]
    d₀, dₖ, dₙ = detectability.(dₘₐₓ, a, [0,collect(1:nₖ-2),nₖ]) 

    ∂S = b - S*(δ+β*sum([I₀, Iₖ...,Iₙ]))
    ∂I₀ = β*S*I₀- (ν + μ + δ + d₀)*I₀
    ∂Iₖ = [β*S*Iₖ[k] - (ν + μ + δ + dₖ[k])*Iₖ[k] + μ*(k > 1 ? Iₖ[k-1] : I₀) for k in 1:nₖ-2]
    ∂Iₙ = β*S*Iₙ - (ν+δ+dₙ)*Iₙ + μ*Iₖ[end]
    ∂R = (ν + d₀)*I₀ + (ν + dₙ)*Iₙ +
        sum([(ν + dₖ[k])*Iₖ[k] for k in 1:nₖ-2]) - δ*R

    return [∂S, ∂I₀, ∂Iₖ..., ∂Iₙ, ∂R]
end


#=
    Parameterize and solve ODE
=#
population_size = 1000.
b, δ, β, ν, μ, a, dₘₐₓ, nₖ = 0.1*population_size, 0.1, 5e-3, 0.2, 0.05, 0.9, 0, 3

# Parameter tuple
θ = (b, δ, β, ν, μ, a, dₘₐₓ, nₖ)

# Initial condition 
u0 = [population_size-1, 1., [0. for i in 2:nₖ]..., 0] 
tspan = (0.,100.)

# Solve ODE 
# (Note: the first time this will take ~10s, and then each subsequent run should be around 3e-3 seconds.
#  this is due to Julia's precompilation).
prob = ODEProblem(SIR, u0, tspan, θ)
@time sol = solve(prob);

#=

    Plot solution 

=#
plt = plot(frame=:box, legend=:outerright)
xlabel!(L"t")
ylabel!("count")
colors = [:darkgreen, :red,  :orange, :yellow, :royalblue]
labs = [L"S", L"I_0", L"I_1", L"I_2", L"R"]
for k in 1:nₖ+2
    plot!(sol.t, [sol.u[t][k] for t in 1:length(sol.u)], lw=2, label=labs[k], c=colors[k])
end
plot!(sol.t, [sum(sol.u[t][2:end-1]) for t in 1:length(sol.u)], lw=2, label=L"\sum_k I_k", c=:black)
plt


#=

    Run and plt heatmap of number of stains (k) and detectability parameter (a)

=#


# This runs on a single core on my laptop in about 10 minutes
function matplt()
    A = 0.0:0.01:1
    K = 3:50
    mat = zeros(length(A), length(K))
    tspan = (0.,1000.)
    popsize = 1000.
    b, δ, β, ν, μ, a, dₘₐₓ= 0.1*popsize, 0.1, 5e-3, 0.2, 0.05, 0.9, 1.0

    @time for (i,a) in enumerate(A)
        @info(a)
        for (j,k) in enumerate(K)
            θ = (b, δ, β, ν, μ, a, dₘₐₓ, k)
            u0 = vcat(popsize-1, 1., zeros(k-1)..., 0) 

            prob = ODEProblem(SIR, u0, tspan, θ)
            sol = solve(prob);
            mat[i,j] = sol.u[end][end]/popsize
        end 
    end
    return mat
end 


mat = matplt()

heatmap(mat, clim=(0.6,.75), size=(800,400),dpi=200,colorbar_title="Equilibrium recovery proportion")
xlabel!(string("Number of variants ", L"(k)"))
ylabel!(string("Detection efficacy decay rate ",L"(a)"))
xticks!(0:5:50, [string(x+3) for x in 0:5:45])

savefig("strains_vs_detectability.png")