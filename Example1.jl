using LinearAlgebra, Distributions, Random, Plots, LaTeXStrings, DataFrames, CSV, Base.Threads
include("SpikeEstimation.jl")
imgFolder = "Figures"
tableFolder = "Tables"

# N = 5000
# σ = sqrt(1.5)
# σ_out = [5,5,4.5]
# σvec = σ^2*ones(N)
# σvec[1:3] = σ_out
# sqrtΣ = Diagonal(sqrt.(σvec))
# k = convert(Int64,floor(log(N)/2))
# jmp = 5
# tol = 3/sqrt(N)
# max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
# vecNbr = 100
# dvec = [0.1,0.5,0.9]
# len_d = length(dvec)
# for m=1:len_d
#     d = dvec[m]
#     M = convert(Int64,ceil(N/d))
#     X = randn(N,M)
#     W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
#     evals = eigvals(W)
#     true_spikes = evals[end:-1:end-2]
#     TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
#     SpikeNbr, SpikeLoc = EstimSpike(TChol,N,c=1)
#     γmin, γplus = EstimSupp(L_list)
#     true_γplus = σ^2*(1+sqrt(d))^2
#     true_γmin = σ^2*(1-sqrt(d))^2
#     x = -0.1+min(γmin,true_γmin):0.001:max(γplus,true_γplus)+0.1
#     yAvrg = EstimDensity(x,L_list,N)
#     MP = x-> true_γmin<x<true_γplus ? sqrt(true_γplus-x)*sqrt(x-true_γmin)/(2*π*d*x*σ^2) : 0
#     true_y = map(x->MP(x),x)
#     p = histogram(evals,bins=evals[1]-0.2:0.1:evals[end]+0.2,normalize=:pdf,label="Empirical ESD")
#     p = plot!(x,yAvrg,linecolor=:red,linewidth=3,label="Estimated Density")
#     p = scatter!(SpikeLoc,0*SpikeLoc,markersize=5,color=:red,marker=:dot,label="Estimated Spikes")
#     p = plot!(x,true_y,linecolor=:blue,linewidth=2,linestyle=:dash,legend=:topright,framestyle=:box,label="True Density")
#     p = scatter!(true_spikes,0*true_spikes,markersize=5,color=:blue,marker=:xcross,label="True Spikes")
#     savefig(p,joinpath(imgFolder, "Ex1Fig"*string(m)*".png"))
#     tb = DataFrame(A=true_spikes,B=SpikeLoc,C=abs.(true_spikes-SpikeLoc))
#     CSV.write(joinpath(tableFolder,"Ex1Tb"*string(m)*".csv"),tb)
# end

# t1 = Threads.@spawn begin
log_file_t1 = open("progress_t1.log", "w")
Nvec = 1000:1000:8000
lenN = length(Nvec)
dvec = [0.1,0.5,0.9]
len_d = length(dvec)
σ = sqrt(1.5)
σ_out = [5,5,4.5]
jmp = 5
vecNbr = 1
SampleNbr = 500 
Percent = zeros(Float64,len_d,lenN)
Avrg = zeros(Float64,len_d,lenN)
@threads for ℓ = 1:len_d
    d = dvec[ℓ]
    @threads for j=1:lenN
        N = Nvec[j]
        M = convert(Int64,ceil(N/d))
        k = convert(Int64,floor(log(N)/2))
        tol = 3/sqrt(N)
        max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
        count = 0
        Sum = 0
        X = zeros(Float64,N,M)
        σvec = σ^2*ones(N)
        σvec[1:3] = σ_out
        sqrtΣ = Diagonal(sqrt.(σvec))
        for ℓ = 1:SampleNbr
            randn!(X)
            W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
            TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
            γmin, γplus = EstimSupp(L_list)
            Nbr, Loc = EstimSpike(TChol,N,c=1.0)
            Sum += Nbr
            if Nbr==3
                count+=1
            end
        end
        msg = "Block 1 with d=$(d) and N=$(N)\n"
        print(msg)
        write(log_file_t1, msg)
        flush(log_file_t1)
        Percent[ℓ,j] = count/SampleNbr
        Avrg[ℓ,j] = Sum/SampleNbr
    end
    tb = DataFrame(A=Nvec,B=Percent[ℓ,:])
    CSV.write(joinpath(tableFolder,"Ex1Prct"*string(ℓ)*".csv"),tb)
    tb = DataFrame(A=Nvec,B=Avrg[ℓ,:])
    CSV.write(joinpath(tableFolder,"Ex1Avrg"*string(ℓ)*".csv"),tb)
end
p4 = plot(Nvec,Percent[1,:],color=:red,linewidth=3,label="",xlabel="N",ylabel="Probability of correct estimation",legend=:topright,framestyle=:box)
p4 = scatter!(Nvec, Percent[1,:], markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
p4 = plot!(Nvec,Percent[2,:],color=:blue,linewidth=3,label="")
p4 = scatter!(Nvec, Percent[2,:], markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
p4 = plot!(Nvec,Percent[3,:],color=:green,linewidth=3,label="")
p4 = scatter!(Nvec, Percent[3,:], markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
savefig(p4,joinpath(imgFolder, "Ex1Fig4.png"))
# end


# t2 = Threads.@spawn begin
# log_file_t2 = open("progress_t2.log", "w")
# Nvec = vcat(200:200:3000,3500:500:8000)
# len_N = length(Nvec)
# dvec = [0.1,0.5,0.9]
# len_d = length(dvec)
# σ = sqrt(1.5)
# σ_out = [5,5,4.5]
# jmp = 5
# vecNbr = 100
# SampleNbr = 200
# SuppErr = zeros(Float64,len_d,SampleNbr,len_N)
# hErr = zeros(Float64,len_d,SampleNbr,len_N)
# @threads for m=1:len_d
#     d = dvec[m]
#     true_γplus = σ^2*(1+sqrt(d))^2
#     true_γmin = σ^2*(1-sqrt(d))^2
#     MP_h = x-> true_γmin<x<true_γplus ? 1/(2*π*d*x*σ^2) : 0
#     @threads for ℓ=1:len_N
#         N = Nvec[ℓ]
#         tol = 3/sqrt(N)
#         k = convert(Int64,floor(log(N)/2))
#         max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
#         M = convert(Int64,ceil(N/d))
#         σvec = σ^2*ones(N)
#         σvec[1:3] = σ_out
#         sqrtΣ = Diagonal(sqrt.(σvec))
#         for κ = 1:SampleNbr
#             X = randn(N,M)
#             W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
#             TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
#             γmin,γplus = EstimSupp(L_list)
#             SuppErr[m,κ,ℓ] = max(abs(true_γmin-γmin),abs(true_γplus-γplus))
#             x = γmin+0.2:0.01:γplus-0.2
#             len_x = length(x)
#             true_h = map(x->MP_h(x),x)
#             yAvrg = EstimDensity(x,L_list,N)
#             hAvrg = zeros(Float64,len_x)
#             for j=1:len_x
#                 hAvrg[j] = yAvrg[j]/(sqrt(γplus-x[j])*sqrt(x[j]-γmin))
#             end
#             hErr[m,κ,ℓ] = maximum(abs.(hAvrg-true_h))
#         end
#         msg = "Block 2 with d=$(d) and N=$(N)\n"
#         print(msg)
#         write(log_file_t2, msg)
#         flush(log_file_t2)
#     end
# end
# close(log_file_t2)
# SuppErr1 = sum(SuppErr[1,:,:],dims=1)/SampleNbr
# SuppErr2 = sum(SuppErr[2,:,:],dims=1)/SampleNbr
# SuppErr3 = sum(SuppErr[3,:,:],dims=1)/SampleNbr
# hErr1 = sum(hErr[1,:,:],dims=1)/SampleNbr
# hErr2 = sum(hErr[2,:,:],dims=1)/SampleNbr
# hErr3 = sum(hErr[3,:,:],dims=1)/SampleNbr
# SuppErr1_std = std(SuppErr[1, :, :], dims=1)
# SuppErr2_std = std(SuppErr[2, :, :], dims=1)
# SuppErr3_std = std(SuppErr[3, :, :], dims=1)
# hErr1_std = std(hErr[1, :, :], dims=1)
# hErr2_std = std(hErr[2, :, :], dims=1)
# hErr3_std = std(hErr[3, :, :], dims=1)
# SuppErr1 = SuppErr1'
# SuppErr2 = SuppErr2'
# SuppErr3 = SuppErr3'
# hErr1 = hErr1'
# hErr2 = hErr2'
# hErr3 = hErr3'
# SuppErr1_std = SuppErr1_std'
# SuppErr2_std = SuppErr2_std'
# SuppErr3_std = SuppErr3_std'
# hErr1_std = hErr1_std'
# hErr2_std = hErr2_std'
# hErr3_std = hErr3_std'
# tb = DataFrame(A=Nvec,B=vec(SuppErr1),C=vec(SuppErr1_std),D=vec(SuppErr2),E=vec(SuppErr2_std),F=vec(SuppErr3),G=vec(SuppErr3_std))
# CSV.write(joinpath(tableFolder,"Ex1SuppErr.csv"),tb)
# tb = DataFrame(A=Nvec,B=vec(hErr1),C=vec(hErr1_std),D=vec(hErr2),E=vec(hErr2_std),F=vec(hErr3),G=vec(hErr3_std))
# CSV.write(joinpath(tableFolder,"Ex1hErr.csv"),tb)
# checkpt = convert(Int64,ceil(len_N/4))
# p2 = plot(Nvec,(SuppErr3[checkpt]*(Nvec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel="Errors in support",legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
# p2 = plot!(Nvec,SuppErr1,yerror=SuppErr1_std,color=:red,linewidth=2,label="")
# p2 = scatter!(Nvec, SuppErr1, markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
# p2 = plot!(Nvec,SuppErr2,yerror=SuppErr2_std,color=:blue,linewidth=2,label="")
# p2 = scatter!(Nvec, SuppErr2, markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
# p2 = plot!(Nvec,SuppErr3,yerror=SuppErr3_std,color=:green,linewidth=2,label="")
# p2 = scatter!(Nvec, SuppErr3, markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
# savefig(p2,joinpath(imgFolder, "Ex1Fig5.png"))
# p3 = plot(Nvec,(hErr1[checkpt]*(Nvec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel=latexstring("\\text{Errors in } \\hat{h}"),legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
# p3 = plot!(Nvec,hErr1,color=:red,linewidth=2,label="")
# p3 = scatter!(Nvec, hErr1, markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
# p3 = plot!(Nvec,hErr2,color=:blue,linewidth=2,label="")
# p3 = scatter!(Nvec, hErr2, markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
# p3 = plot!(Nvec,hErr3,color=:green,linewidth=2,label="")
# p3 = scatter!(Nvec, hErr3, markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
# savefig(p3,joinpath(imgFolder, "Ex1Fig6.png"))
# end

# t3 = Threads.@spawn begin
# Nvec = 100:100:8000 
# len_N = length(Nvec)
# d = 0.5
# σ = sqrt(1.5)
# σ_out = [5,5,4.5]
# jmp = 5
# vecNbr = 1
# EigTime = zeros(Float64,len_N)
# Time = zeros(Float64,len_N)
# for ℓ=1:len_N
#     N = Nvec[ℓ]
#     tol = 3/sqrt(N)
#     k = convert(Int64,floor(log(N)/2))
#     max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
#     M = convert(Int64,ceil(N/d))
#     X = randn(N,M)
#     σvec = σ^2*ones(N)
#     σvec[1:3] = σ_out
#     sqrtΣ = Diagonal(sqrt.(σvec))
#     W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
#     tmp = @timed begin
#         eigvals(W) 
#     end
#     EigTime[ℓ] = tmp.time
#     tmp = @timed begin
#         TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
#         Nbr,Loc = EstimSpike(TChol,N,c=1)
#     end
#     Time[ℓ] += tmp.time
# end
# p4 = plot(Nvec,EigTime,color=:red,linewidth=2,label="",xlabel="N",ylabel="Time (in seconds)",legend=:topleft,framestyle=:box)
# p4 = scatter!(Nvec, EigTime, markersize=4, color=:red, marker=:diamond, label="Eigenvalue Computation")
# p4 = plot!(Nvec,Time,color=:orange,linewidth=2,label="")
# p4 = scatter!(Nvec, Time, markersize=4, color=:orange, marker=:square, label="Lanczos Approach")
# savefig(p4,joinpath(imgFolder, "Ex1Fig7.png"))
# end

# fetch(t1)
# fetch(t2)
# fetch(t3)