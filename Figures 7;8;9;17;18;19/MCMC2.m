%%% Author Sergios Agapiou

%%% Code corresponding to Compartmental Model 1 in Section 2.4.1
%%% of "Modeling of Covid-19 Pandemic in Cyprus" by Agapiou et al.
%%% Code based on the code of "Substantial undocumented infection
%%% facilitates the rapid dissemination of novel coronavirus (sars-cov-2)",
%%% in Science 368, by Li et al.
%%% This code is simplified by not considering the metapopulation
%%% structure. An independence sampler is used instead of Kalman filter
%%% techniques.

%%% We consider 6 periods of 14 days, and estimate alpha and beta (and
%%% hence R) in each of them. 

%%% This script considers all incidents (as opposed to local transmission
%%% only)


%%% Function MCMC2 takes as input:
%%% alpha0_vec: previous alpha chain
%%% initial_state: state at the end of the previous period
%%% initial_incidence: incidence for all days, even in the future, as
%%% determined by the end of the last period
%%% hypalp1, hypalp2: hyperparameters for Beta prior on alpha
%%% hypbet1, hypbet2: hyperparameters for Gamma prior on beta
%%% init_day: initial day of period
%%% final_day: final day of period
%%% old_iday: initial day of previous period
%%% period: index of current period.

%%% Function MCMC2 runs an independence sampler and returns chains of
%%% alpha, beta, incidence, states and R

function [state_vec, alpha_vec, beta_vec,incidence_vec, Re_vec]=MCMC2(alpha0_vec,initial_state,initial_incidence,hypalp1, hypalp2, hypbet1, hypbet2, init_day, final_day, old_iday, period)

%%% Other parameter values needed in the model 
mu=0.5;     % relative transmissibility
theta=0;    % movement factor
Z=5.1;      % latency period
D=3.5;      % infectious period
pop=875000; % population of Cyprus

%%% Number of days in period under consideration
num_days=final_day-init_day+1;

%%% Initial states and incidence
S0=initial_state(1,:);
E0=initial_state(2,:);
Is0=initial_state(3,:);
Ia0=initial_state(4,:);
y0=initial_incidence;

%%% Number of iterations of the independence sampler
iter=length(initial_state(1,:));

%%% Define vector for alpha, beta and final state chains
alpha_vec=zeros(1,iter);
beta_vec=zeros(1,iter);
state_vec=zeros(4,iter); %to store final states

%%% Load incidence data
incidence_real= readtable('cyprus_data_official.csv'); 
incidence_real = incidence_real.Incidence;
incidence_real = incidence_real(1:end); 
incidence_real1 = incidence_real(init_day:final_day); % incidence data for this period



%%% Likelihood variance as in Li et al 2020
OEV=max(4,(incidence_real1.^2)./(4)); 

%%% Initialize Markov chains by drawing from respective priors and
%%% computing resulting state and incidence vectors
alpha_vec(1)=betarnd(hypalp1,hypalp2);
beta_vec(1)=gamrnd(hypbet1,hypbet2); 

%%% After the first period, the initial numbers of susceptibles, exposed, 
%%% infected-reported, infected-unreported, are drawn by sampling the final
%%% states at end of the previous period, subject to a small correction 
%%% to close any accumulated gaps when comparing total real incidence and 
%%% total sampled incidence at the end of the last period.  
u0=randi([1,iter]);
sgap=sum(incidence_real(1:(init_day-1)))-sum(y0(1:(init_day-1),u0));
sa=alpha0_vec(u0); %a=Ir/(Ir+Iu) aIr+aIu=Ir Iu=(1-a)Ir/a I=Ir/a
gap=incidence_real(init_day-1)-y0(init_day-1,u0);
a=alpha0_vec(u0);  %a=Ir/(Ir+Iu) aIr+aIu=Ir Iu=(1-a)Ir/a I=Ir/a
[incidence_vec(:,1),state_vec(:,1)]=predict_cyprus_MCMC_mult(num_days,initial_incidence(:,u0),alpha_vec(1),beta_vec(1),S0(u0)-sgap/sa,E0(u0)+gap/a,Is0(u0)+gap,Ia0(u0)+gap*(1-a)/a,0,init_day);
nsamp=0;

%%% Independence Sampler
tic
for i=1:iter-1
    alpha=betarnd(hypalp1,hypalp2);%betarnd(hypalp1,hypalp2);%0.5+rand/2;%betarnd(hypalp1,hypalp2);
    beta=gamrnd(hypbet1,hypbet2);
    u0=randi([1,iter]);
    sgap=sum(incidence_real(1:(init_day-1)))-sum(y0(1:(init_day-1),u0));
    sa=alpha0_vec(u0);
    gap=incidence_real(init_day-1)-y0(init_day-1,u0);
    a=alpha0_vec(u0);    %Ir/(Ir+Iu)=a, hence I=Ir+Iu=Ir/a and Iu=Ir/a-Ir=Ir(1/a-1)
    [incidence_prop, state_prop]=predict_cyprus_MCMC_mult(num_days,initial_incidence(:,u0),alpha,beta,S0(u0)-sgap/sa,E0(u0)+gap/a,Is0(u0)+gap,Ia0(u0)+gap*(1-a)/a,0,init_day);
    aux=(sum(((incidence_vec(init_day:(init_day+num_days-1),i)-incidence_real1).^2)./OEV)-sum(((incidence_prop(init_day:(init_day+num_days-1))-incidence_real1).^2)./OEV))/2;
    u=rand;
    while (log(u)>aux)
        alpha=betarnd(hypalp1,hypalp2);%betarnd(hypalp1,hypalp2);%0.5+rand/2;%betarnd(hypalp1,hypalp2);
        beta=gamrnd(hypbet1,hypbet2);
        u0=randi([1,iter]);
        sgap=sum(incidence_real(1:(init_day-1)))-sum(y0(1:(init_day-1),u0));
        sa=alpha0_vec(u0);
        gap=incidence_real(init_day-1)-y0(init_day-1,u0);
        a=alpha0_vec(u0);   
        [incidence_prop, state_prop]=predict_cyprus_MCMC_mult(num_days,initial_incidence(:,u0),alpha,beta,S0(u0)-sgap/sa,E0(u0)+gap/a,Is0(u0)+gap,Ia0(u0)+gap*(1-a)/a,0,init_day);
        aux=(sum(((incidence_vec(init_day:(init_day+num_days-1),i)-incidence_real1).^2)./OEV)-sum(((incidence_prop(init_day:(init_day+num_days-1))-incidence_real1).^2)./OEV))/2;
        u=rand;
        
        nsamp=nsamp+1;
    end
    alpha_vec(i+1)=alpha;
    beta_vec(i+1)=beta;
    incidence_vec(:,i+1)=incidence_prop;
    state_vec(:,i+1)=state_prop;   
    i
end
toc
nsamp

%%% Compute R chain from alpha and beta chains
Re_vec=alpha_vec.*beta_vec*D + (1-alpha_vec).*mu.*beta_vec*D;


















