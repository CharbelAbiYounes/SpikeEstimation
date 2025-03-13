using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Distributed")

f = open("hosts.txt")
nodes = readlines(f)
close(f)

using LinearAlgebra, Distributions, Random, Plots, LaTeXStrings, DataFrames, CSV, Distributed
imgFolder = "Figures"
tableFolder = "Tables"

num_procs = 10
addprocs([nodes[2] for j in 1:num_procs],tunnel=true)
addprocs([nodes[3] for j in 1:num_procs],tunnel=true)
addprocs([nodes[4] for j in 1:num_procs],tunnel=true)
addprocs(num_procs-1)

@everywhere begin
    using LinearAlgebra, Distributions, Random
    include("SpikeEstimation.jl")
end

log_file_t2 = open("progress_t2.log", "w")
log_file_t2Det = open("progress_t2Det.log", "w")

Nvec = vcat(200:200:3000,3500:500:6000)
len_N = length(Nvec)
# dvec = [0.1,0.5,0.9]
dvec = [0.5,0.9]
len_d = length(dvec)
σ = sqrt(1.5)
σ_out = [5,5,4.5]
jmp = 5
vecNbr = 100
SampleNbr = 200

SuppErr = zeros(Float64,len_d,SampleNbr,len_N)
hErr = zeros(Float64,len_d,SampleNbr,len_N)

@everywhere function process_sample(d, N, σ, σ_out, jmp, vecNbr)
    σvec = σ^2 * ones(N)
    σvec[1:3] = σ_out
    sqrtΣ = Diagonal(sqrt.(σvec))

    tol = 3 / sqrt(N)
    k = convert(Int64, floor(log(N) / 2))
    max_iter = convert(Int64, ceil(max(6 * log(N) + 24, sqrt(N))))
    M = convert(Int64, ceil(N / d))

    X = randn(N, M)
    W = sqrtΣ * (1 / M * X * X') * sqrtΣ' |> Symmetric
    TChol, L_list = CholList(W, tol, k, jmp, max_iter, vecNbr)
    γmin, γplus = EstimSupp(L_list)

    true_γplus = σ^2 * (1 + sqrt(d))^2
    true_γmin = σ^2 * (1 - sqrt(d))^2
    supp_err = max(abs(true_γmin - γmin), abs(true_γplus - γplus))

    x = γmin + 0.2:0.01:γplus - 0.2
    len_x = length(x)
    MP_h = x -> true_γmin < x < true_γplus ? 1 / (2 * π * d * x * σ^2) : 0
    true_h = map(MP_h, x)

    yAvrg = EstimDensity(x, L_list, N)
    hAvrg = yAvrg ./ (sqrt.(γplus .- x) .* sqrt.(x .- γmin))
    h_err = maximum(abs.(hAvrg - true_h))

    return (supp_err, h_err)
end

for m in 1:len_d
    d = dvec[m]
    for ℓ in 1:len_N
        N = Nvec[ℓ]

        results = pmap(_ -> process_sample(d, N, σ, σ_out, jmp, vecNbr), 1:SampleNbr;batch_size=10)

        SuppErr[m, :, ℓ] .= first.(results)
        hErr[m, :, ℓ] .= last.(results)

        msg = "Block 2 with d=$(d) and N=$(N)\n"
        print(msg)
        write(log_file_t2, msg)
        flush(log_file_t2)
    end
end

close(log_file_t2)
close(log_file_t2Det)

SuppErr1 = sum(SuppErr[1,:,:],dims=1)/SampleNbr
SuppErr2 = sum(SuppErr[2,:,:],dims=1)/SampleNbr
# SuppErr3 = sum(SuppErr[3,:,:],dims=1)/SampleNbr
hErr1 = sum(hErr[1,:,:],dims=1)/SampleNbr
hErr2 = sum(hErr[2,:,:],dims=1)/SampleNbr
# hErr3 = sum(hErr[3,:,:],dims=1)/SampleNbr
SuppErr1_std = std(SuppErr[1, :, :], dims=1)
SuppErr2_std = std(SuppErr[2, :, :], dims=1)
# SuppErr3_std = std(SuppErr[3, :, :], dims=1)
hErr1_std = std(hErr[1, :, :], dims=1)
hErr2_std = std(hErr[2, :, :], dims=1)
# hErr3_std = std(hErr[3, :, :], dims=1)
SuppErr1 = SuppErr1'
SuppErr2 = SuppErr2'
# SuppErr3 = SuppErr3'
hErr1 = hErr1'
hErr2 = hErr2'
# hErr3 = hErr3'
SuppErr1_std = SuppErr1_std'
SuppErr2_std = SuppErr2_std'
# SuppErr3_std = SuppErr3_std'
hErr1_std = hErr1_std'
hErr2_std = hErr2_std'
# hErr3_std = hErr3_std'
# tb = DataFrame(A=Nvec,B=vec(SuppErr1),C=vec(SuppErr1_std),D=vec(SuppErr2),E=vec(SuppErr2_std),F=vec(SuppErr3),G=vec(SuppErr3_std))
# CSV.write(joinpath(tableFolder,"Ex1SuppErr.csv"),tb)
tb = DataFrame(A=Nvec,B=vec(SuppErr1),C=vec(SuppErr1_std),D=vec(SuppErr2),E=vec(SuppErr2_std))
CSV.write(joinpath(tableFolder,"Ex1SuppErr0509S.csv"),tb)
# tb = DataFrame(A=Nvec,B=vec(hErr1),C=vec(hErr1_std),D=vec(hErr2),E=vec(hErr2_std),F=vec(hErr3),G=vec(hErr3_std))
# CSV.write(joinpath(tableFolder,"Ex1hErr.csv"),tb)
tb = DataFrame(A=Nvec,B=vec(hErr1),C=vec(hErr1_std),D=vec(hErr2),E=vec(hErr2_std))
CSV.write(joinpath(tableFolder,"Ex1hErrd0509S.csv"),tb)
checkpt = convert(Int64,ceil(len_N/4))
# p2 = plot(Nvec,(SuppErr3[checkpt]*(Nvec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel="Errors in support",legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
p2 = plot(Nvec,(SuppErr2[checkpt]*(Nvec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel="Errors in support",legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
p2 = plot!(Nvec,SuppErr1,yerror=SuppErr1_std,color=:red,linewidth=2,label="")
p2 = scatter!(Nvec, SuppErr1, markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
p2 = plot!(Nvec,SuppErr2,yerror=SuppErr2_std,color=:blue,linewidth=2,label="")
p2 = scatter!(Nvec, SuppErr2, markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
# p2 = plot!(Nvec,SuppErr3,yerror=SuppErr3_std,color=:green,linewidth=2,label="")
# p2 = scatter!(Nvec, SuppErr3, markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
savefig(p2,joinpath(imgFolder, "Ex1Fig5d0509S.png"))
# p3 = plot(Nvec,(hErr1[checkpt]*(Nvec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel=latexstring("\\text{Errors in } \\hat{h}"),legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
p3 = plot(Nvec,(hErr1[checkpt]*(Nvec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel=latexstring("\\text{Errors in } \\hat{h}"),legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
p3 = plot!(Nvec,hErr1,color=:red,linewidth=2,label="")
p3 = scatter!(Nvec, hErr1, markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
p3 = plot!(Nvec,hErr2,color=:blue,linewidth=2,label="")
p3 = scatter!(Nvec, hErr2, markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
# p3 = plot!(Nvec,hErr3,color=:green,linewidth=2,label="")
# p3 = scatter!(Nvec, hErr3, markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
savefig(p3,joinpath(imgFolder, "Ex1Fig6d0509S.png"))