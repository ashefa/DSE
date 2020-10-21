include("PowerFlowRes57.jl")
include("case57.m")
E=PF[:,2].*cos((π/180)PF[:,3])+PF[:,2].*sin((π/180)PF[:,3])im
Iₗ=complex(zeros(size(branch,1),2))
for i in 1:size(branch,1)
    Iₗ[i,1]=((PF1[i,4]/baseMVA)+(PF1[i,5]/baseMVA)im)/E[Int(PF1[i,2])]
    Iₗ[i,2]=((PF1[i,6]/baseMVA)+(PF1[i,7]/baseMVA)im)/E[Int(PF1[i,3])]
end

PMU=[1, 4, 6, 9, 15, 20, 24, 28, 30, 32, 36, 38, 41, 46, 50, 53, 57]
# PMU=[2, 6, 7, 9] #case14
# PMU=[1, 2, 3, 4, 8, 10, 11, 15, 16, 17, 21, 23] #case24_ieee_rts
# PMU=[2, 5, 9, 12, 15, 17, 21, 25, 29, 34, 37, 42, 45, 49, 53, 56,
# 62, 63, 68, 70, 71, 75, 77, 80, 85, 86, 91, 94, 102, 105, 110, 114] #case118
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
z=B*E+(2/3000)*B*E.*(rand(size(B,1))+rand(size(B,1))im)
M=(transpose(B)*(W^(-1))*B)^(-1)*transpose(B)*(W^(-1))
x̂=M*z