%%% Code integrating the Compartmental Model 1 in Section 2.4.1
%%% of "Modeling of Covid-19 Pandemic in Cyprus" by Agapiou et al.
%%% Code based on the code of "Substantial undocumented infection
%%% facilitates the rapid dissemination of novel coronavirus (sars-cov-2)",
%%% in Science 368, by Li et al.
%%% This code is simplified by not considering the metapopulation
%%% structure. 


function x=SEIR(x)

dt=1;
tmstep=1;
pop=875000;
%integrate forward for one day

%S,E,Is,Ia,obs,...,beta,mu,theta,Z,alpha,D
Sidx=1;
Eidx=2;
Isidx=3;
Iaidx=4;
obsidx=5;
betaidx=6;
muidx=7;
thetaidx=8;
Zidx=9;
Didx=10;
alphaidx=11;

S=zeros(1,tmstep+1);
E=zeros(1,tmstep+1);
Is=zeros(1,tmstep+1);
Ia=zeros(1,tmstep+1);
Incidence=zeros(1,tmstep+1);
obs=0;

%initialize S,E,Is,Ia and parameters
S(1)=x(1);
E(1)=x(2);
Is(1)=x(3);
Ia(1)=x(4);
beta=x(6);
mu=x(7);
theta=x(8);
Z=x(9);
D=x(10);
alpha=x(11);

%start integration
tcnt=0;
for t=dt:dt:tmstep
    tcnt=tcnt+1;
    dt1=dt;
    
    %first step
    Eexps=dt1*(beta)*S(tcnt)*Is(tcnt)/pop;
    Eexpa=dt1*(mu)*(beta)*S(tcnt)*Ia(tcnt)/pop;
    Einfs=dt1*(alpha).*E(tcnt)./(Z);
    Einfa=dt1*((1-alpha)).*E(tcnt)./(Z);
    Erecs=dt1*Is(tcnt)./(D);
    Ereca=dt1*Ia(tcnt)./(D);
    
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);% making sure no negatives
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    
    %%%%%%%%%%stochastic version
    Eexps=poissrnd(Eexps); % U1
    Eexpa=poissrnd(Eexpa); % U2
    Einfs=poissrnd(Einfs); % U5
    Einfa=poissrnd(Einfa); % U6
    Erecs=poissrnd(Erecs); % U9
    Ereca=poissrnd(Ereca); % U10

    sk1=-Eexps-Eexpa; 
    ek1=Eexps+Eexpa-Einfs-Einfa;
    isk1=Einfs-Erecs;
    iak1=Einfa-Ereca;
    ik1i=Einfs; % This is the new reported incidents (without removing the recovered ones)
    
    
    %second step
    Ts1=S(tcnt)+sk1/2;
    Te1=E(tcnt)+ek1/2;
    Tis1=Is(tcnt)+isk1/2;
    Tia1=Ia(tcnt)+iak1/2;
        
    Eexps=dt1*(beta)*Ts1*Tis1/pop;
    Eexpa=dt1*(mu)*(beta)*Ts1*Tia1/pop;
    Einfs=dt1*(alpha)*Te1/(Z);
    Einfa=dt1*((1-alpha))*Te1/(Z);
    Erecs=dt1*Tis1/(D);
    Ereca=dt1*Tia1/(D);
    
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);% making sure no negatives
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    
    %%%%%%%%%%stochastic version
    Eexps=poissrnd(Eexps);
    Eexpa=poissrnd(Eexpa);
    Einfs=poissrnd(Einfs);
    Einfa=poissrnd(Einfa);
    Erecs=poissrnd(Erecs);
    Ereca=poissrnd(Ereca);

    sk2=-Eexps-Eexpa;
    ek2=Eexps+Eexpa-Einfs-Einfa;
    isk2=Einfs-Erecs;
    iak2=Einfa-Ereca;
    ik2i=Einfs;
    
    %third step
    Ts2=S(tcnt)+sk2/2;
    Te2=E(tcnt)+ek2/2;
    Tis2=Is(tcnt)+isk2/2;
    Tia2=Ia(tcnt)+iak2/2;
    
     
    Eexps=dt1*(beta)*Ts2*Tis2/pop;
    Eexpa=dt1*(mu)*(beta)*Ts2*Tia2/pop;
    Einfs=dt1*(alpha)*Te2/(Z);
    Einfa=dt1*((1-alpha))*Te2/(Z);
    Erecs=dt1*Tis2/(D);
    Ereca=dt1*Tia2/(D);
    
  
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);% making sure no negatives
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    
    %%%%%%%%%%stochastic version
    Eexps=poissrnd(Eexps);
    Eexpa=poissrnd(Eexpa);
    Einfs=poissrnd(Einfs);
    Einfa=poissrnd(Einfa);
    Erecs=poissrnd(Erecs);
    Ereca=poissrnd(Ereca);

    sk3=-Eexps-Eexpa;
    ek3=Eexps+Eexpa-Einfs-Einfa;
    isk3=Einfs-Erecs;
    iak3=Einfa-Ereca;
    ik3i=Einfs;
    
    %fourth step
    Ts3=S(tcnt)+sk3;
    Te3=E(tcnt)+ek3;
    Tis3=Is(tcnt)+isk3;
    Tia3=Ia(tcnt)+iak3;
    
    
    Eexps=dt1*(beta)*Ts3*Tis3/pop;
    Eexpa=dt1*(mu)*(beta)*Ts3*Tia3/pop;
    Einfs=dt1*(alpha)*Te3/(Z);
    Einfa=dt1*((1-alpha))*Te3/(Z);
    Erecs=dt1*Tis3/(D);
    Ereca=dt1*Tia3/(D);
    
    
    Eexps=max(Eexps,0);Eexpa=max(Eexpa,0);% making sure no negatives
    Einfs=max(Einfs,0);Einfa=max(Einfa,0);
    Erecs=max(Erecs,0);Ereca=max(Ereca,0);
    
    %%%%%%%%%%stochastic version
  
    Eexps=poissrnd(Eexps);
    Eexpa=poissrnd(Eexpa);
    Einfs=poissrnd(Einfs);
    Einfa=poissrnd(Einfa);
    Erecs=poissrnd(Erecs);
    Ereca=poissrnd(Ereca);

    sk4=-Eexps-Eexpa;
    ek4=Eexps+Eexpa-Einfs-Einfa;
    isk4=Einfs-Erecs;
    iak4=Einfa-Ereca;
    ik4i=Einfs;
    
    %%%%%
    S(tcnt+1)=S(tcnt)+round(sk1/6+sk2/3+sk3/3+sk4/6);
    E(tcnt+1)=E(tcnt)+round(ek1/6+ek2/3+ek3/3+ek4/6);
    Is(tcnt+1)=Is(tcnt)+round(isk1/6+isk2/3+isk3/3+isk4/6);
    Ia(tcnt+1)=Ia(tcnt)+round(iak1/6+iak2/3+iak3/3+iak4/6);
    Incidence(tcnt+1)=round(ik1i/6+ik2i/3+ik3i/3+ik4i/6);
    obs=Incidence(tcnt+1);
end
%%%update x
x(Sidx)=S(tcnt+1);
x(Eidx)=E(tcnt+1);
x(Isidx)=Is(tcnt+1);
x(Iaidx)=Ia(tcnt+1);
x(obsidx)=obs;
