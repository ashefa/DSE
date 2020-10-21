include("PowerFlowRes14.jl")
include("case14.m")
E=PF[:,2].*cos((π/180)PF[:,3])+PF[:,2].*sin((π/180)PF[:,3])im
Iₗ=complex(zeros(size(branch,1),2))
for i in 1:size(branch,1)
    Iₗ[i,1]=((PF1[i,4]/baseMVA)+(PF1[i,5]/baseMVA)im)/E[Int(PF1[i,2])]
    Iₗ[i,2]=((PF1[i,6]/baseMVA)+(PF1[i,7]/baseMVA)im)/E[Int(PF1[i,3])]
end

PMU=[2, 6, 7, 9]
r=branch[:,3];x=branch[:,4];y=1./(r+x*im);b=branch[:,5];y₀=b/2*im;
nLinePMU=0
for k in 1:size(branch,1)
    for k′ in 1:size(PMU,1)
        if branch[k,1]==PMU[k′]||branch[k,2]==PMU[k′]
            for k′′ in 1:size(PMU,1)
                if branch[k,1]==PMU[k′]&&branch[k,2]==PMU[k′′]
                    nLinePMU=nLinePMU+1
                end
            end
        end
    end
end
A=zeros(size(branch,1)+nLinePMU,size(bus,1))
Y=complex(zeros(size(branch,1)+nLinePMU,size(branch,1)+nLinePMU));
yₛ=complex(zeros(size(branch,1)+nLinePMU,size(bus,1)));y₀=[y₀;zeros(nLinePMU,1)]
I₎=complex(zeros(size(branch,1)+nLinePMU,1));
a′=0
for k in 1:size(branch,1)
    for k′ in 1:size(PMU,1)
        if branch[k,1]==PMU[k′]||branch[k,2]==PMU[k′]
            if branch[k,1]==PMU[k′]&&branch[k,2]!=PMU[k′]
                A[k+a′,PMU[k′]]=1
                A[k+a′,Int(branch[k,2])]=-1
                Y[k+a′,k+a′]=y[k]
                yₛ[k+a′,PMU[k′]]=y₀[k]
                I₎[k+a′]=Iₗ[k,1]
            elseif branch[k,2]==PMU[k′]&&branch[k,1]!=PMU[k′]
                A[k+a′,PMU[k′]]=1
                A[k+a′,Int(branch[k,1])]=-1
                Y[k+a′,k+a′]=y[k]
                yₛ[k+a′,PMU[k′]]=y₀[k]
                I₎[k+a′]=Iₗ[k,2]
            end
            for k′′ in 1:size(PMU,1)
                if branch[k,1]==PMU[k′]&&branch[k,2]==PMU[k′′]
                    A[k+a′,PMU[k′]]=1
                    A[k+a′,Int(branch[k,2])]=-1
                    Y[k+a′,k+a′]=y[k]
                    yₛ[k+a′,PMU[k′]]=y₀[k]
                    I₎[k+a′]=Iₗ[k,1]
                    A[k+a′+1,PMU[k′′]]=1
                    A[k+a′+1,Int(branch[k,1])]=-1
                    Y[k+a′+1,k+a′+1]=y[k]
                    yₛ[k+a′+1,PMU[k′′]]=y₀[k]
                    I₎[k+a′+1]=Iₗ[k,2]
                    a′=a′+1
                end
            end
        end
    end
end
A′=A
A=A[find(A[i,:]!=zeros(size(A,2)) for i in 1:size(A,1)),:]
Y=Y[find(Y[i,:]!=zeros(size(Y,2)) for i in 1:size(Y,1)),find(Y[i,:]!=zeros(size(Y,2)) for i in 1:size(Y,1))]
yₛ=yₛ[find(A′[i,:]!=zeros(size(A′,2)) for i in 1:size(A′,1)),find(A′[:,i]!=zeros(size(A′,1)) for i in 1:size(A′,2))]
ΙΙ=zeros(size(PMU,1),size(bus,1))
a′′=0
for i in 1:size(bus,1)
    for j in 1:size(PMU,1)
        if i==PMU[j]
            ΙΙ[1+a′′,i]=1
            a′′=a′′+1
        end
    end
end
B=[ΙΙ;Y*A+yₛ]
W=1000*eye(size(B,1))
NS = size(bus,1); T = 24; i = 1:NS; t = 1:T;
z=B*repmat(E,1,T)+(2/3000)*B*repmat(E,1,T).*(rand(size(B,1),T)+rand(size(B,1),T)im)
M=(transpose(B)*(W^(-1))*B)^(-1)*transpose(B)*(W^(-1))
#x̂=M*z
W_p=1000*eye(size(B,2))
K=[W zeros(size(W,1),size(W_p,1));zeros(size(W_p,1),size(W,1)) W_p]
B_p=[B;eye(size(B,2))]
M_p=(transpose(B_p)*(K^(-1))*B_p)^(-1)*transpose(B_p)*(K^(-1))
NS = size(bus,1); T = 24; i = 1:NS; t = 1:T;
x̃=complex(zeros(NS,T))
Sₚ=complex(zeros(NS,T))
bₚ=complex(zeros(NS,T))
α=repmat([0.8],1,NS);β=repmat([0.5],1,NS)
z_p=[z;x̃];x̄=repmat(E,1,T)+(2/3000)*repmat(E,1,T).*(rand(size(E,1),T)+rand(size(E,1),T)im)
x̂=complex(zeros(NS,T))
x̂₁=complex(zeros(NS,T))

for t in 1:T
    if t in 1:3
        Sₚ[i,t]=x̄[i,t]
        bₚ[i,t]=x̄[i,t+1]-x̄[i,t]
        x̃[i,t]=x̄[i,t]
    else
        Sₚ[i,t-1]=α[i].*(x̂[i,t-1])+(1-α[i]).*x̃[i,t-1]
        bₚ[i,t-1]=β[i].*(Sₚ[i,t-1]-Sₚ[i,t-2])+(1-β[i]).*(bₚ[i,t-2])
        x̃[i,t]=α[i].*(1+β[i]).*(x̂[i,t-1])+
        (1+β[i]).*(1-α[i]).*x̃[i,t-1]-β[i].*Sₚ[i,t-2]+(1-β[i]).*bₚ[i,t-2]
    end
    z_p[:,t]=[z[:,t];x̃[:,t]] 
    x̂[:,t]=M_p*z_p[:,t]
    x̂₁[:,t]=M*z[:,t]
end

MSE=sum((real(x̂[i,:])-repmat(real(E),1,T)[i,:]).^2+(imag(x̂[i,:])-repmat(imag(E),1,T)[i,:]).^2 for i in 1:NS)/2NS
MSE₁=sum((real(x̂₁[i,:])-repmat(real(E),1,T)[i,:]).^2+(imag(x̂₁[i,:])-repmat(imag(E),1,T)[i,:]).^2 for i in 1:NS)/2NS

    using PyPlot
# plot line 1
plot(1:24,MSE ,label="DSE")
# plot line 2
plot(1:24,MSE₁, label="SSE")
# set x and y axis labels
xlabel("T"); ylabel(L"$MSE$")
legend() # turn on legend
# add a title
title(Comparison Between SSE & DSE)