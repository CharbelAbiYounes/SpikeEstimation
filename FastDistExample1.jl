using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Distributed")

using LinearAlgebra, Distributions, Random, Plots, LaTeXStrings, DataFrames, CSV, Distributed
imgFolder = "Figures"
tableFolder = "Tables"

f = open("hosts.txt")
nodes = readlines(f)
close(f)

num_procs = 2
addprocs([nodes[2] for j in 1:num_procs],tunnel=true)
addprocs([nodes[3] for j in 1:num_procs],tunnel=true)
addprocs([nodes[4] for j in 1:num_procs],tunnel=true)
addprocs(num_procs-1)

@everywhere begin
    using LinearAlgebra, Distributions, Random
    include("SpikeEstimation.jl")
end

logfile = open("progress_Distt1.log", "w")
Detlogfile = open("progress_Distt1Det.log", "w")

# Nvec = 100:100:8000
Nvec = 7900:8000
lenN = length(Nvec)
dvec = [0.1,0.5,0.9]
len_d = length(dvec)
σ = sqrt(1.5)
σout = [5,5,4.5]
jmp = 5
vecNbr = 1
SampleNbr = 200 
Percent = zeros(Float64,len_d,lenN)
SpikeNbr = zeros(Float64,len_d,lenN)
LanTime = zeros(Float64,len_d,lenN)
EigTime = zeros(Float64,len_d,lenN)

@everywhere function process_sample(N, M, sqrtΣ, tol, k, jmp, max_iter, vecNbr)
    X = randn(N,M)
    W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
    EigDuration = @timed begin
        evals = eigvals(W)
    end
    LanDuration = @timed begin
        TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
        Nbr, Loc = EstimSpike(TChol,N,c=1.0)
    end
    return (EigDuration.time, LanDuration.time, Nbr)
end

for m in 1:len_d
    d = dvec[m]
    for ℓ in 1:lenN
        N = Nvec[ℓ]
        M = convert(Int64,ceil(N/d))
        k = convert(Int64,floor(log(N)/2))
        tol = 3/sqrt(N)
        max_iter = convert(Int64,ceil(max(6*log(N)+24,sqrt(N))))
        σvec = σ^2*ones(N)
        σvec[1:3] = σout
        sqrtΣ = Diagonal(sqrt.(σvec))
        results = pmap(_-> process_sample(N, M, sqrtΣ, tol, k, jmp, max_iter, vecNbr), 1:SampleNbr)

        EigTime[m, ℓ] = sum(first.(results)) / SampleNbr
        LanTime[m, ℓ] = sum(x -> x[2], results) / SampleNbr  
        SpikeNbr[m, ℓ] = sum(x -> x[3], results) / SampleNbr
        Percent[m, ℓ] = sum(Int(x[3] == 3) for x in results) / SampleNbr  
        msg = "d=$(d) and N=$(N)\n"
        print(msg)
        write(logfile, msg)
        flush(logfile)
        tb = DataFrame(A=Nvec[1:ℓ],B=Percent[m,1:ℓ])
        CSV.write(joinpath(tableFolder,"L010509Prct"*string(m)*".csv"),tb)
        tb = DataFrame(A=Nvec[1:ℓ],B=SpikeNbr[m,1:ℓ])
        CSV.write(joinpath(tableFolder,"L010509Avrg"*string(m)*".csv"),tb)
        tb = DataFrame(A=Nvec[1:ℓ],B=LanTime[m,1:ℓ])
        CSV.write(joinpath(tableFolder,"L010509LanTime"*string(m)*".csv"),tb)
        tb = DataFrame(A=Nvec[1:ℓ],B=EigTime[m,1:ℓ])
        CSV.write(joinpath(tableFolder,"L010509EigTime"*string(m)*".csv"),tb)
    end
end
close(logfile)
close(Detlogfile)

# p = plot(Nvec,Percent[1,:],color=:red,linewidth=3,label="",xlabel="N",ylabel="Probability of correct estimation",legend=:topright,framestyle=:box)
# p = scatter!(Nvec, Percent[1,:], markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
# p = plot!(Nvec,Percent[2,:],color=:blue,linewidth=3,label="")
# p = scatter!(Nvec, Percent[2,:], markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
# p = plot!(Nvec,Percent[3,:],color=:green,linewidth=3,label="")
# p = scatter!(Nvec, Percent[3,:], markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
# savefig(p,joinpath(imgFolder, "Avrg.png"))
# p = plot(Nvec,EigTime[1,:],color=:red,linewidth=2,label="",xlabel="N",ylabel="Time (in seconds)",legend=:topleft,framestyle=:box)
# p = scatter!(Nvec, EigTime[1,:], markersize=4, color=:red, marker=:diamond, label="Eigenvalue Computation")
# p = plot!(Nvec,LanTime[1,:],color=:orange,linewidth=2,label="")
# p = scatter!(Nvec, LanTime[1,:], markersize=4, color=:orange, marker=:square, label="Lanczos Approach")
# savefig(p,joinpath(imgFolder, "Time01.png"))
# p = plot(Nvec,EigTime[2,:],color=:red,linewidth=2,label="",xlabel="N",ylabel="Time (in seconds)",legend=:topleft,framestyle=:box)
# p = scatter!(Nvec, EigTime[2,:], markersize=4, color=:red, marker=:diamond, label="Eigenvalue Computation")
# p = plot!(Nvec,LanTime[2,:],color=:orange,linewidth=2,label="")
# p = scatter!(Nvec, LanTime[2,:], markersize=4, color=:orange, marker=:square, label="Lanczos Approach")
# savefig(p,joinpath(imgFolder, "Time05.png"))
# p = plot(Nvec,EigTime[3,:],color=:red,linewidth=2,label="",xlabel="N",ylabel="Time (in seconds)",legend=:topleft,framestyle=:box)
# p = scatter!(Nvec, EigTime[3,:], markersize=4, color=:red, marker=:diamond, label="Eigenvalue Computation")
# p = plot!(Nvec,LanTime[3,:],color=:orange,linewidth=2,label="")
# p = scatter!(Nvec, LanTime[3,:], markersize=4, color=:orange, marker=:square, label="Lanczos Approach")
# savefig(p,joinpath(imgFolder, "Time05.png"))