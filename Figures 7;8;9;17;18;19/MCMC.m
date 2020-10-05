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

clear all;
close all;
format long;


%%% Priors on reporting rate alpha and transmission rate beta are Beta and
%%% Gamma respectively. Hyperparameters for the 2 priors for the 1st period are as below.
hypalp1=2;
hypalp2=2;
hypbet1=3/2;
hypbet2=3/2;

%%% Outbreak initialized by drawing initial uniform integers between 0 and maxE0=maxIa0=10 
%%% for the number of exposed and the number of infected-unreported on 4
%%% Mar 2020, that is 3 days before the first recorded incident.
maxE0=10;
maxIa0=10;

%%% Other parameter values needed in the model 
mu=0.5;     % relative transmissibility
theta=0;    % movement factor
Z=5.1;      % latency period
D=3.5;      % infectious period
pop=875000; % population of Cyprus


%%% Number of iterations of the Independence Sampler (for each of the 6
%%% time periods)
iter=100; 


%%% Load incidence data
incidence_real= readtable('cyprus_data_official.csv');
incidence_real = incidence_real.Incidence;
incidence_real = incidence_real(4:end-5); % incidence data for this period


num_days=length(incidence_real); % total number of days in data
num_times=14; % 14 days in each period under examination 
num_periods=6; % 6 fortnight periods under examination
incidence_real1 = incidence_real(1:num_times); % incidence data for 1st fornight period


%%% Define Markov chain vectors, for first fortnight period
alpha_vec=zeros(num_periods,iter);   % reporting rate alpha chain
beta_vec=zeros(num_periods,iter);    % tansimission rate beta chain 
Re_vec=zeros(num_periods,iter);      % effective reproduction number R chain
Ia0_vec=zeros(1,iter);               % initial number of unreported infections chain
E0_vec=zeros(1,iter);                % initial number of exposed chain

%%% Define vectors for incidence and states
incidence_vec=zeros(num_days,iter,num_periods);        % to store timeseries of incidences in each period
state_vec=zeros(4,iter,num_periods);                   % to store final states
init_incidence=zeros(1,num_days);                      % to store initial incidence
%%% Initial incidence is used to implement the delay mechanism in recording
%%% of an infected-reported case, see predict_cyprus_MCMC.m
mean_incidence_vec=zeros(num_times,num_periods);       % to store mean incidence
credible_intervals_vec=zeros(2,num_times,num_periods); % to store credible intervals for incidence


%%% Likelihood variance as in Li et al 2020
OEV1=max(4,(incidence_real1.^2)./4);


%%% Initialize Markov chains by drawing from respective priors and
%%% computing resulting state and incidence vectors
alpha_vec(1,1)=betarnd(hypalp1,hypalp2);   %initialize alpha chain
beta_vec(1,1)=gamrnd(hypbet1,hypbet2);     %initialize beta chain
E0_vec(1)=randi([0,maxE0]);                %initialize E0 chain
Ia0_vec(1)=randi([0,maxIa0]);              %initialize Ia0 chain
[incidence_vec(:,1,1),state_vec(:,1,1)]=predict_cyprus_MCMC(num_times,init_incidence,alpha_vec(1,1)',beta_vec(1,1)',pop-E0_vec(1)-Ia0_vec(1),E0_vec(1),0,Ia0_vec(1),0);


%%% Independence Sampler 
nsamp=0;
tic
for i=1:iter-1
    alpha=betarnd(hypalp1,hypalp2);
    beta=gamrnd(hypbet1,hypbet2);
    E0=randi([0,maxE0]);
    Ia0=randi([0,maxIa0]);
    [incidence_prop, state_prop]=predict_cyprus_MCMC(num_times,init_incidence,alpha,beta,pop-E0-Ia0,E0,0,Ia0,0);
    aux=(sum(((incidence_vec(1:num_times,i)-incidence_real1).^2)./OEV1)-sum(((incidence_prop(1:num_times)'-incidence_real1).^2)./OEV1))/2;
    u=rand;
    while (log(u)>aux)
        alpha=betarnd(hypalp1,hypalp2);
        beta=gamrnd(hypbet1,hypbet2);
        E0=randi([0,maxE0]);
        Ia0=randi([0,maxIa0]);
        [incidence_prop, state_prop]=predict_cyprus_MCMC(num_times,init_incidence,alpha,beta,pop-E0-Ia0,E0,0,Ia0,0);
        aux=(sum(((incidence_vec(1:num_times,i)-incidence_real1).^2)./OEV1)-sum(((incidence_prop(1:num_times)'-incidence_real1).^2)./OEV1))/2;
        u=rand;
        nsamp=nsamp+1;
    end
    alpha_vec(1,i+1)=alpha;
    beta_vec(1,i+1)=beta;
    Ia0_vec(i+1)=Ia0;
    E0_vec(i+1)=E0;
    incidence_vec(:,i+1,1)=incidence_prop;
    state_vec(:,i+1,1)=state_prop; 
    i
end
toc
nsamp

%%% Compute R chain from alpha and beta chains
Re_vec(1,:)=alpha_vec(1,:).*beta_vec(1,:)*D + (1-alpha_vec(1,:)).*mu.*beta_vec(1,:)*D;


%%% Feed resulting data as initial values for running independence samplers
%%% for later fortnight periods. Independence sampler implemented in
%%% MCMC2.m. 
%%% The Beta prior on alpha, is chosen increasingly skewed towards higher
%%% reporting rates.
for j=2:num_periods
        [state_vec(:,:,j), alpha_vec(j,:), beta_vec(j,:),incidence_vec(:,:,j), Re_vec(j,:)]=MCMC2(alpha_vec(j-1,:),state_vec(:,:,j-1),incidence_vec(:,:,j-1),max(2,(j+1)/2),2,1,1,(j-1)*num_times+1,(j)*num_times,(j-2)*num_times+1,j);
end

