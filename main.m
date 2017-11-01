% This code is based on the follwing paper
% Inversion-Based Output Tracking and Unknown Input Reconstruction of Square Discrete-Time Linear Systems By E.Naderi and K.Khorasani
% This is an illustrative code. It is neither optimized nor generalized
% This code does not cover all aspects of the paper, for instance the case system has zero on the unit circle 
% This code is solely developed for my own research
% Author: Esmaeil Naderi

clc
clear

%% Example 1 of the paper
Ts=1;
sys=zpk([1.5 0.5],[0 0],1,Ts); 
sys=ss(sys);
%sys=canon(sys,'companion');

%% Example 2 of the paper
% A=[0.6 -0.3 0  0; 0.1 1.0  0  0;-0.4  -1.5 0.4  -0.3;...
%     0.3   1.1   0.2  0.9];
% B=[ 0  0.4;0   0;0 -0.1;0.1   0.1];
% 
% C=[1 2 3 4;2 1 5 6];
% sys=ss(A,B,C,0,Ts);

%% Initialization
n=size(sys.a,1);
l=size(sys.c,1);
m=size(sys.b,2);

ST=300; %% simulation time
nD=20;  %% delay

%% Get the system output

x0=rand(n,1);
t=0:Ts:ST*Ts;
u=randn(m,numel(t));u=u';
[y,t,x]=lsim(sys,u,t,x0);x=x'; % Note this is x(k) not x^{(1)}(k)=Q'*x(k)
ug=[t u];
yg=[t y];

%%  Calculate M, Ahat and F


Dn=compute_HL(sys,n-1);
Cn=compute_OL(sys,n-1);
In=[eye(m) zeros(m,n*m-m)];
MZ=sys.a-sys.b*In*pinv(Dn)*Cn;
[Vs,Ds]=eig(MZ');

M0=[];
Ahat=[];
k=1;

while k <=n
    if abs(Ds(k,k))<1
        if imag(Ds(k,k))<1e-2
            Ds(k,k)=real(Ds(k,k));
            Vs(:,k)=real(Vs(:,k));
            M0=[M0 Vs(:,k)];
            As=size(Ahat,1);
            Ahat=[Ahat zeros(As,1);zeros(1,As) Ds(k,k)];
        else
            M0=[M0 (Vs(:,k)+Vs(:,k+1))/2 (Vs(:,k)-Vs(:,k+1))/(2i)];
            As=size(Ahat,1);
            Ahat=[Ahat zeros(As,2);zeros(2,As) [real(Ds(k,k)) ...
                imag(Ds(k,k));-imag(Ds(k,k)) real(Ds(k,k))]];
            k=k+1;
        end
    end
    k=k+1;
end
M0=M0';
Ahat=Ahat';
F=M0*sys.b*In*pinv(Dn);
%%  Uncomment this part for the case that M0 is rank deficient
% Uncomment this part if rank rank(M0)-(n-\beta)=1 and the first row of 
% the M0 is linearly dependent to other rows. Coding for a general case
% according to the paper is beyond the scope of this illusterative code.

% Mn=Mn(2:end,:);
% Ahat=Ahat(2:end,2:end);
% F=F(2:end,:);
% 
% rM=1;
% rA=0.1*rand;
% DC=(eye(n*l)-Dn*pinv(Dn))*Cn;
% Kn=rand(rM,n*l);
% KDC=Kn*DC;
% Mp=sylvester(rA,-MZ,-KDC);
% Fp=Mp*sys.b*In*pinv(Dn)+Kn*(eye(n*m)-Dn*pinv(Dn));
% M0=[Mp; M0];
% Ahat=[rA 0 0;zeros(2,1) Ahat];
% F=[Fp ;F];

%%  Calculate the parameters of the UIO and FIR filter
% Note that this code should be supported by many if then conditions to 
% be in accordance of the paper, but this is just an illustration. 

mr=rank(M0);

[Q,R]=qr(M0'); L=R'; % LQ Decomposition
Abar=Q'*sys.a*Q; 
Bbar=Q'*sys.b;
Cbar=sys.c*Q;
sysMn=ss(Abar,Bbar,Cbar,sys.d,Ts); % transforming system with T=Q'

% Partitioning system matrices
A11=sysMn.a(1:mr,1:mr);
A12=sysMn.a(1:mr,mr+1:end);
A21=sysMn.a(mr+1:end,1:mr);
A22=sysMn.a(mr+1:end,mr+1:end);
B1=sysMn.b(1:mr,:);
B2=sysMn.b(mr+1:end,:);
C1=sysMn.c(:,1:mr);
C2=sysMn.c(:,mr+1:end);

% Defining Mq and inv(Mq)
Mq=L(:,1:mr);
iMq=pinv(Mq);

% We left the condition that B1 is not full column rank. You can code it
% using the material of the paper.

B1s=size(B1,2);
B2s=size(B2,2);
if rank(B1)==B1s 
    Az=A22-B2*pinv(B1)*A12;
    iAz=inv(Az);
    Bz=[B2*pinv(B1) A21-B2*pinv(B1)*A11];
    Cz1=[sysMn.d*pinv(B1) sysMn.c(:,1:mr)-sysMn.d*pinv(B1)*A11];
    Cz2=sysMn.c(:,mr+1:end)-sysMn.d*pinv(B1)*A12;
end
     

%% Integration for calculating unknown states and inputs

unknownInput=[];
MPstate=[];
NMPstate=[];
Theta=[];

zeta0=rand(mr,1); % Initialize the UIO
xnmp=rand(n-mr,1); %Initialize FIR filter

for j=1:ST   
    if j>n
        Y=Build_Mat(yg,j-n,1,n); % build Yn
    else 
        Y=zeros(n*l,1);
    end
    
    zeta1=Ahat*zeta0+F*Y; %integrate UIO
    x11=iMq*zeta1;
    x10=iMq*zeta0;
    MPstate=[MPstate x11];

    
    beta=Bz*[x11;x10];    
    Theta=[Theta beta];
    if j>nD+n %integrate FIR filter
        for h=1:nD-1
            xnmp=iAz*(xnmp-Theta(:,j-h+1));
            uk=pinv(B1)*(MPstate(:,j-h+1)-A11*MPstate(:,j-h)-A12*xnmp);
        end
    else
        uk=zeros(m,1);
    end
    unknownInput=[unknownInput uk];
    NMPstate=[NMPstate xnmp];

    zeta0=zeta1;
end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This plot is only working for this example.
%%%  You have to change it for your example
         
xq=Q'*x;
plot(xq(1,1:50))
hold on
plot(MPstate(2:51),'o')
ylabel('x_1')
legend('system state','Reconstructed state')
title('Minimum Phase States Reconstruction')

figure

plot(xq(2,1:50))
hold on
plot(-nD+2:1:70-nD+1, NMPstate(2:71),'o')
ylabel('x_1')
legend('system state','Reconstructed state')
title('Non-Minimum Phase States Reconstruction (Note the shift)')

figure
plot(u(2:50))
hold on
plot(-nD+1:1:70-nD, unknownInput(2:71),'o')
legend('system input','Reconstructed input')
title('Unknown Input Reconstruction (Note the shift)')
         
         
%% End of code

