using LinearAlgebra, Distributions, Plots, LaTeXStrings, DataFrames, CSV
include("SpikeEstimation")
imgFolder = "Figures"
tableFolder = "Tables"

N = 5000
d = 0.1
M = convert(Int64,ceil(N/d))
X = randn(N,M)
σ = sqrt(1.5)
σ_out = [4,4,3.5]
σvec = σ^2*ones(N)
σvec[1:3] = σ_out
sqrtΣ = Diagonal(sqrt.(σvec))
W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
evals = eigvals(W)
true_spikes = evals[end:-1:end-2]
k = convert(Int64,floor(log(N)/2))
jmp = 5
tol = 2.5/sqrt(N)
max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
vecNbr = 100
TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
SpikeNbr, SpikeLoc = EstimSpike(TChol,N)
γmin, γplus = EstimSupp(L_list)
x = -0.1+γmin:0.001:γplus+0.1
yAvrg = EstimDensity(x,L_list,N)
true_γplus = σ^2*(1+sqrt(d))^2
true_γmin = σ^2*(1-sqrt(d))^2
MP = x-> true_γmin<x<true_γplus ? sqrt(true_γplus-x)*sqrt(x-true_γmin)/(2*π*d*x*σ^2) : 0
true_y = map(x->MP(x),x)
p1 = histogram(evals,bins=γmin-0.2:0.1:γplus+0.2,normalize=:pdf,label="Empirical ESD")
p1 = plot!(x,yAvrg,linecolor=:red,linewidth=3,label="Estimated Density")
p1 = scatter!(SpikeLoc,0*SpikeLoc,markersize=5,color=:red,marker=:dot,label="Estimated Spikes")
p1 = plot!(x,true_y,linecolor=:blue,linewidth=2,linestyle=:dash,legend=:topright,framestyle=:box,label="True Density")
p1 = scatter!(true_spikes,0*true_spikes,markersize=5,color=:blue,marker=:xcross,label="True Spikes")
savefig(p1,joinpath(imgFolder, "Ex1Fig1.png"))
tb1 = DataFrame(true_spikes,SpikeLoc,abs.(true_spikes-SpikeLoc))
CSV.write(joinpath(tableFolder,"Ex1Tb1.csv"),tb1)

d = 0.5
M = convert(Int64,ceil(N/d))
X = randn(N,M)
σvec = σ^2*ones(N)
σvec[1:3] = σ_out
sqrtΣ = Diagonal(sqrt.(σvec))
W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
evals = eigvals(W)
true_spikes = evals[end:-1:end-2]
TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
SpikeNbr, SpikeLoc = EstimSpike(TChol,N)
γmin, γplus = EstimSupp(L_list)
x = -0.1+γmin:0.001:γplus+0.1
yAvrg = EstimDensity(x,L_list,N)
true_γplus = σ^2*(1+sqrt(d))^2
true_γmin = σ^2*(1-sqrt(d))^2
MP = x-> true_γmin<x<true_γplus ? sqrt(true_γplus-x)*sqrt(x-true_γmin)/(2*π*d*x*σ^2) : 0
true_y = map(x->MP(x),x)
p2 = histogram(evals,bins=γmin-0.2:0.1:γplus+0.2,normalize=:pdf,label="Empirical ESD")
p2 = plot!(x,yAvrg,linecolor=:red,linewidth=3,label="Estimated Density")
p2 = scatter!(SpikeLoc,0*SpikeLoc,markersize=5,color=:red,marker=:dot,label="Estimated Spikes")
p2 = plot!(x,true_y,linecolor=:blue,linewidth=2,linestyle=:dash,legend=:topright,framestyle=:box,label="True Density")
p2 = scatter!(true_spikes,0*true_spikes,markersize=5,color=:blue,marker=:xcross,label="True Spikes")
savefig(p2,joinpath(imgFolder, "Ex1Fig2.png"))
tb2 = DataFrame(true_spikes,SpikeLoc,abs.(true_spikes-SpikeLoc))
CSV.write(joinpath(tableFolder,"Ex1Tb2.csv"),tb2)

d = 0.9
M = convert(Int64,ceil(N/d))
X = randn(N,M)
σvec = σ^2*ones(N)
σvec[1:3] = σ_out
sqrtΣ = Diagonal(sqrt.(σvec))
W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
evals = eigvals(W)
true_spikes = evals[end:-1:end-2]
k = convert(Int64,floor(log(N)/2))
TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
SpikeNbr, SpikeLoc = EstimSpike(TChol,N)
γmin, γplus = EstimSupp(L_list)
x = -0.1+γmin:0.001:γplus+0.1
yAvrg = EstimDensity(x,L_list,N)
true_γplus = σ^2*(1+sqrt(d))^2
true_γmin = σ^2*(1-sqrt(d))^2
MP = x-> true_γmin<x<true_γplus ? sqrt(true_γplus-x)*sqrt(x-true_γmin)/(2*π*d*x*σ^2) : 0
true_y = map(x->MP(x),x)
p3 = histogram(evals,bins=γmin-0.2:0.1:γplus+0.2,normalize=:pdf,label="Empirical ESD")
p3 = plot!(x,yAvrg,linecolor=:red,linewidth=3,label="Estimated Density")
p3 = scatter!(SpikeLoc,0*SpikeLoc,markersize=5,color=:red,marker=:dot,label="Estimated Spikes")
p3 = plot!(x,true_y,linecolor=:blue,linewidth=2,linestyle=:dash,legend=:topright,framestyle=:box,label="True Density")
p3 = scatter!(true_spikes,0*true_spikes,markersize=5,color=:blue,marker=:xcross,label="True Spikes")
savefig(p3,joinpath(imgFolder, "Ex1Fig3.png"))
tb3 = DataFrame(true_spikes,SpikeLoc,abs.(true_spikes-SpikeLoc))
CSV.write(joinpath(tableFolder,"Ex1Tb3.csv"),tb3)

Nvec = 100:100:8000
lenN = length(Nvec)
dvec = [0.1,0.5,0.9]
len_d = length(dvec)
σ = sqrt(1.5)
σ_out = [4,4,3.5]
jmp = 5
vecNbr = 1
SampleNbr = 200 
Percent = zeros(Float64,len_d,lenN)
for ℓ = 1:len_d
    d = dvec[ℓ]
    for j=1:lenN
        N = Nvec[j]
        M = convert(Int64,ceil(N/d))
        k = convert(Int64,floor(log(N)/2))
        tol = 2.5/sqrt(N)
        max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
        count = 0
        for ℓ = 1:SampleNbr
            X = randn(N,M)
            σvec = σ^2*ones(N)
            σvec[1:3] = σ_out
            sqrtΣ = Diagonal(sqrt.(σvec))
            W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
            TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
            γmin, γplus = EstimSupp(L_list)
            Nbr, Loc = EstimSpike(TChol,N)
            if Nbr==3
                count+=1
            end
        end
        Percent[ℓ,j] = count/SampleNbr
    end
end
p4 = plot(Nvec,Percent[1,:],color=:red,linewidth=3,label="",xlabel="N",ylabel="Probability of correct estimation",legend=:topright,framestyle=:box)
p4 = scatter!(Nvec, Percent[1,:], markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
p4 = plot(Nvec,Percent[2,:],color=:blue,linewidth=3,label="")
p4 = scatter!(Nvec, Percent[2,:], markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
p4 = plot(Nvec,Percent[3,:],color=:green,linewidth=3,label="")
p4 = scatter!(Nvec, Percent[3,:], markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
savefig(p4,joinpath(imgFolder, "Ex1Fig4.png"))

Nvec = 100:100:8000 
len_N = length(Nvec)
dvec = [0.1,0.5,0.9]
len_d = length(dvec)
σ = sqrt(1.5)
σ_out = [4,4,3.5]
jmp = 5
vecNbr = 100
SampleNbr = 200
SuppErr = zeros(Float64,len_d,SampleNbr,len_N)
hErr = zeros(Float64,len_d,SampleNbr,len_N)
for m=1:len_d
    d = dvec[m]
    true_γplus = σ^2*(1+sqrt(d))^2
    true_γmin = σ^2*(1-sqrt(d))^2
    x = true_γmin+0.25:0.01:true_γplus-0.25
    len_x = length(x)
    MP_h = x-> true_γmin<x<true_γplus ? 1/(2*π*d*x*σ^2) : 0
    true_h = map(x->MP_h(x),x)
    for ℓ=1:len_N
        N = Nvec[ℓ]
        tol = 2.5/sqrt(N)
        k = convert(Int64,floor(log(N)/2))
        max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
        M = convert(Int64,ceil(N/d))
        σvec = σ^2*ones(N)
        σvec[1:3] = σ_out
        sqrtΣ = Diagonal(sqrt.(σvec))
        for κ = 1:SampleNbr
            X = randn(N,M)
            W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
            TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
            γmin,γplus = EstimSupp(L_list)
            yAvrg = EstimDensity(x,L_list,N)
            hAvrg = zeros(Float64,len_x)
            for j=1:len_x
                hAvrg[j] = yAvrg[j]/(sqrt(γplus-x[j])*sqrt(x[j]-γmin))
            end
            SuppErr[m,κ,ℓ] = max(abs(true_γmin-γmin),abs(true_γplus-γplus))
            hErr[m,κ,ℓ] = maximum(abs.(hAvrg-true_h))
        end
    end
end
SuppErr1 = sum(SuppErr[1,:,:],dims=1)/SampleNbr
SuppErr2 = sum(SuppErr[2,:,:],dims=1)/SampleNbr
SuppErr3 = sum(SuppErr[3,:,:],dims=1)/SampleNbr
hErr1 = sum(hErr[1,:,:],dims=1)/SampleNbr
hErr2 = sum(hErr[2,:,:],dims=1)/SampleNbr
hErr3 = sum(hErr[3,:,:],dims=1)/SampleNbr
SuppErr1 = SuppErr1'
SuppErr2 = SuppErr2'
SuppErr3 = SuppErr3'
hErr1 = hErr1'
hErr2 = hErr2'
hErr3 = hErr3'
checkpt = convert(Int64,ceil(len_N/4))
p5 = plot(Nvec,(SuppErr1[checkpt]*(N_vec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel="Errors in support",legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
p5 = plot!(Nvec,SuppErr1,color=:red,linewidth=2,label="")
p5 = scatter!(Nvec, SuppErr1, markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
p5 = plot!(Nvec,SuppErr2,color=:blue,linewidth=2,label="")
p5 = scatter!(Nvec, SuppErr2, markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
p5 = plot!(Nvec,SuppErr3,color=:green,linewidth=2,label="")
p5 = scatter!(Nvec, SuppErr3, markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
savefig(p5,joinpath(imgFolder, "Ex1Fig5.png"))
p6 = plot(Nvec,(hErr1[checkpt]*(N_vec[checkpt])^(1/2))*(Nvec).^(-1/2),color=:orange,linewidth=3,label=L"\mathrm{O}(N^{-1/2})",xlabel="N",ylabel="Errors in ̂"*L"\hat{h}",legend=:topright,yscale=:log10,linestyle=:dash,framestyle=:box)
p6 = plot!(Nvec,(hErr1[checkpt]*(N_vec[checkpt])^(1/6))*(Nvec).^(-1/6),color=:purple,linewidth=3,label=L"\mathrm{O}(N^{-1/6})")
p6 = plot!(Nvec,hErr1,color=:red,linewidth=2,label="")
p6 = scatter!(Nvec, hErr1, markersize=4, color=:red, marker=:diamond, label="d="*string(dvec[1]))
p6 = plot!(Nvec,hErr2,color=:blue,linewidth=2,label="")
p6 = scatter!(Nvec, hErr2, markersize=4, color=:blue, marker=:square, label="d="*string(dvec[2]))
p6 = plot!(Nvec,hErr3,color=:green,linewidth=2,label="")
p6 = scatter!(Nvec, hErr3, markersize=4, color=:green, marker=:circ, label="d="*string(dvec[3]))
savefig(p6,joinpath(imgFolder, "Ex1Fig6.png"))

Nvec = 100:100:8000 
len_N = length(Nvec)
d = 0.9
σ = sqrt(1.5)
σ_out = [4,4,3.5]
jmp = 5
vecNbr = 1
EigTime = zeros(Float64,len_N)
Time = zeros(Float64,len_N)
for ℓ=1:len_N
    N = Nvec[ℓ]
    tol = 2.5/sqrt(N)
    k = convert(Int64,floor(log(N)/2))
    max_iter = convert(Int64,ceil(max(6*log(N)+24,N/4,sqrt(N))))
    M = convert(Int64,ceil(N/d))
    X = randn(N,M)
    σvec = σ^2*ones(N)
    σvec[1:3] = σ_out
    sqrtΣ = Diagonal(sqrt.(σvec))
    W = sqrtΣ*(1/M*X*X')*sqrtΣ'|>Symmetric
    tmp = @timed begin
        eigvals(W) 
    end
    EigTime[ℓ] = tmp.time
    tmp = @timed begin
        TChol,L_list = CholList(W,tol,k,jmp,max_iter,vecNbr)
        Nbr,Loc = EstimSpike(TChol,N)
    end
    Time[ℓ] += tmp.time
end
p7 = plot(Nvec,EigTime,color=:red,linewidth=2,label="",xlabel="N",ylabel="Time (in seconds)",legend=:topleft,framestyle=:box)
p7 = scatter!(Nvec, EigTime, markersize=4, color=:red, marker=:diamond, label="Eigenvalue Computation")
p7 = plot!(Nvec,Time,color=:orange,linewidth=2,label="")
p7 = scatter!(Nvec, Time, markersize=4, color=:orange, marker=:square, label="Lanczos Approach")
savefig(p7,joinpath(imgFolder, "Ex1Fig7.png"))