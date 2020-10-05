%%% Code corresponding to Compartmental Model 1 in Section 2.4.1
%%% of "Modeling of Covid-19 Pandemic in Cyprus" by Agapiou et al.
%%% Code based on the code of "Substantial undocumented infection
%%% facilitates the rapid dissemination of novel coronavirus (sars-cov-2)",
%%% in Science 368, by Li et al.

%%% This function requires initial values S0, E0, IS0, IA0, y0, num_times and
%%% init_incidence. 
%%% The function returns synthetic incedence data for each of the days in
%%% the period under consideration produced in the following way. 
%%% We evolve the state day by day. For each day, for each
%%% infected-reported case, we draw a random delay and record this case in
%%% a suitable later period. Initial incidence is used to carry forward
%%% cases that although have appeared as infected in previous periods, they
%%% need to be recorded during this period, due to the random delay.
%%% The function also returns the state at the end of the period under
%%% consideration.

function [incidence_synthetic,state]=predict_cyprus_MCMC(num_times,init_incidence,alpha,beta,S0,E0,Is0,Ia0,y0)



pop = 875000; %load population
Td=6;         %average reporting delay
a=1.85;       %shape parameter of gamma distribution
b=Td/a;       %scale parameter of gamma distribution

%observation operator: obs=Hx
H=zeros(1,5+6);
H(1,5)=1; 

mu=0.5;   %relative transmissibility
theta=0;  %movement factor
Z=5.1;    %latency period
D=3.5;    %infectious period

total_days=length(init_incidence);


x=[S0,E0,Is0,Ia0,y0,beta,mu,theta,Z,D,alpha]'; %initialize x

x_mat=zeros(size(x,1),num_times); %matrix in which we store x for all times
obs_temp=init_incidence;          %records of observations

for t=1:num_times
    x=SEIR_cyprus(x);
    obs_cnt=H*x;%new infection
    %add reporting delay
    if obs_cnt>0
       rnd=ceil(gamrnd(a,b,obs_cnt,1));
       for h=1:length(rnd)
           if (t+rnd(h)<=total_days)
               obs_temp(t+rnd(h))=obs_temp(t+rnd(h))+1;
           end
       end
    end
    x_mat(:,t) = x;
end

state=x_mat(1:4,num_times);
incidence_synthetic = obs_temp; 


