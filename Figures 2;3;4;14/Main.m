clear;close all;clc;
 %----------------------------
Npop=875000; % popoulation of Cyprus
load data_cyprus 
%------------------------------------
% the following if different data set is used
%------------------------------------
% Confirmed_cy=data_cyprus(:,1);
% Recovered_cy=data_cyprus(:,2);
% Deaths_cy=data_cyprus(:,3);
%----------------------------------
time=datetime(datetime('1-March-2020'):1:datetime('31-May-2020'));
% %------------------------------------
% % Data preparation
% %------------------------------------
Confirmed=cumsum(Confirmed_cy);
Recovered=cumsum(Recovered_cy);
Deaths=cumsum(Deaths_cy);
minNum= min(20,round(0.25*max(Confirmed))); 
if isempty(Confirmed_cy)
    warning('"Confirmed" is an empty array. Check the value of "minNum". Computation aborted.')
    return
end

Recovered(Confirmed<=minNum)=[]; %Remove cases where only few infectious are recorded (to avoid bad
%    initial conditions)
Deaths(Confirmed<=minNum)=[];
time(Confirmed<=minNum)= [];
Confirmed(Confirmed<=minNum)=[];
%----------------------------------
% % plot of Confirmed vs Recovered
% %----------------------------------
% figure (1)
% clf
% semilogy(time,Confirmed_cy,'bo',time,Recovered_cy+0.1,'ro');
% ylabel('Number of cases')
% xlabel('time (days)')
% leg = {'Confirmed (reported)', 'Recovered (reported)'};
% legend(leg{:},'location','southoutside')
% set(gcf,'color','w')
% grid on
% axis tight
% set(gca,'yscale','lin')
% %-------------------------------------
% % Active case
% %--------------------------
% figure (2)
% clf
% semilogy(time,Confirmed_cy-Recovered_cy-Deaths_cy,'bo');
% ylabel('Number of cases')
% xlabel('time (days)')
% leg = { 'Active (reported)'};
% legend(leg{:},'location','southoutside')
% set(gcf,'color','w')
% grid on
% axis tight
% set(gca,'yscale','lin')
%--------------------------------------
% Initial values
%--------------------------------------
E0 = Confirmed(1); % Initial number of exposed cases. Unknown but unlikely to be zero.
I0 = Confirmed(1); % Initial number of infectious cases. Unknown but unlikely to be zero.
Q0 = Confirmed(1)-Recovered(1)-Deaths(1);
R0 = Recovered(1);
D0 = Deaths(1);
%%------------------
%%This part studies the influence of γ to the rest of parameters 
%%------------------
lambda_guess = [0.01,0.001,10]; % recovery rate
kappa_guess = [0.01,0.001,1]; % death rate
alpha_guess=0.7
beta_guess=0.2;
Q_guess=0.5;
Alpha=[];
Beta=[];
Gamma=[];
Delta=[];

LT_guess=[1:1:7];
for ii=1:length(LT_guess)
lT_guess= LT_guess(ii) % latent time in days
guess = [alpha_guess,beta_guess,1/lT_guess,Q_guess,lambda_guess,kappa_guess];
Active = Confirmed-Recovered-Deaths;
Active(Active<0) = 0; % No negative number possible Parameter estimation with the lsqcurvefit function
[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,lambdaFun,kappaFun] = ...
    fit_SEIQRDP_fix(Active,Recovered,Deaths,Npop,E0,I0,time,guess);%,'Display','off');
Alpha=[Alpha;alpha1];
Beta=[Beta;beta1];
Gamma=[Gamma;gamma1];
Delta=[Delta;delta1];

dt = 0.1; % time step
time1 = datetime(time(1)):dt:datetime(datestr(floor(datenum('31-May-2020'))+datenum(0))); % predicts for up to 0 days
N = numel(time1);
t = [0:N-1].*dt;

 % Call of the function SEIQRDP.m with the fitted parameters
[S,E,I,Q,R,D,P] = SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,E0,I0,Q0,R0,D0,t,lambdaFun,kappaFun);

% figure(4)
% semilogy(time1(601:end),Q(601:end));
% hold on
% semilogy(time,Active,'bo');
% leg = {' Estimated','Observed' };
% legend(leg{:},'location','northeast')
% ylabel('Active cases')
% xlabel('time (days)')
% set(gcf,'color','w')
% grid on
% axis tight
% set(gca,'yscale','lin')

figure(5)
semilogy(time1(601:end),R(601:end));
hold on
semilogy(time,Recovered,'bo');
leg = {'Estimated','Observed' };
legend(leg{:},'location','northeast')
ylabel('Recovered cases')
xlabel('time (days)')
set(gcf,'color','w')
grid on
axis tight
set(gca,'yscale','lin')

figure(6)
semilogy(time1,E+I);
hold on
ylabel('Exposed + Infected cases')
xlabel('time (days)')
leg = {'Estimated' };
legend(leg{:},'location','northeast')
set(gcf,'color','w')
grid on
axis tight
set(gca,'yscale','lin')
end

figure(3) % to produce these figures need to go to fit_SEIQRDP_fix which keeps the value of γ fixed
clf
semilogy(1./Gamma,Alpha,'bo',1./Gamma,Beta,'r<',1./Gamma, Delta,'k>');
ylabel('')
xlabel('Latent time \gamma^{-1}(days)')
leg = {'\zeta', '\beta', '\delta^{-1}'};
legend(leg{:},'location','northeast')
set(gcf,'color','w')
grid on
axis tight
set(gca,'yscale','lin')
axis([0 7 0 1])
%---------------------
% This part tries to choose the right value of γ 
%------------------
lambda_guess = [0.01,0.001,10]; % recovery rate
kappa_guess = [0.01,0.001,1]; % death rate
alpha_guess=0.7
beta_guess=0.2;
Q_guess=0.5;
Alpha=[];
Beta=[];
Gamma=[];
Delta=[];

LT_guess=[1:1:7];
for ii=1:length(LT_guess)
lT_guess= LT_guess(ii) % latent time in days
guess = [alpha_guess,beta_guess,1/lT_guess,Q_guess,lambda_guess,kappa_guess];
Active = Confirmed-Recovered-Deaths;
Active(Active<0) = 0; % No negative number possible Parameter estimation with the lsqcurvefit function
[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,lambdaFun,kappaFun] = ...
    fit_SEIQRDP(Active,Recovered,Deaths,Npop,E0,I0,time,guess);%,'Display','off');
Alpha=[Alpha;alpha1];
Beta=[Beta;beta1];
Gamma=[Gamma;gamma1];
Delta=[Delta;delta1];
end
%-----------------------------------------------
%
% Fit and prediction part starts here
%
%-----------------------------------------------
lambda_guess = [0.01,0.001,10]; % recovery rate
kappa_guess = [0.01,0.001,1]; % death rate
alpha_guess=0.8;
beta_guess=0.3;
Q_guess=0.5;
lT_guess= 1/3; % latent time in days^(-1)
guess = [alpha_guess,beta_guess,1/lT_guess,Q_guess,lambda_guess,kappa_guess];
Active = Confirmed-Recovered-Deaths;
Active(Active<0) = 0; % No negative number possible
% Parameter estimation with the lsqcurvefit function
Alpha_t=[];
Beta_t=[];
Gamma_t=[];
Delta_t=[];
RE=[];
RE_R=[];
RE_D=[];
dates_pr=[0 7 16 44 59 ]; % which correspond to 31/5 24/5 15/5, 17/4, 2/4 
for iii=1:length(dates_pr)
    Nest=dates_pr(iii)
[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,lambdaFun,kappaFun] = ...
    fit_SEIQRDP(Active(1:length(Active)-Nest),Recovered(1:length(Active)-Nest),...
    Deaths(1:length(Active)-Nest),Npop,E0,I0,time(1:length(Active)-Nest),guess);%,'Display','off');
Alpha_t=[Alpha_t;alpha1];
Beta_t=[Beta_t;beta1];
Gamma_t=[Gamma_t;gamma1];
Delta_t=[Delta_t;delta1];

%-------------------------------------------------------------
%
% Simulate the epidemic outbreak based on the fitted parameters
%
%--------------------------------------------------------------
dt = 0.1; % time step
time1 = datetime(time(1)):dt:datetime(datestr(floor(datenum('31-May-2020'))+datenum(30))); % predicts for up to 30 days
N = numel(time1);
t = [0:N-1].*dt;

% recovery and death rate 

% lambda_f=Lambda1(1)./(1+exp(-Lambda1(2)*(t-Lambda1(3))));
% kappa_f=Kappa1(1).*exp(-Kappa1(2).*(t-Kappa1(3)).^2);
% figure(1)
% plot(t,lambda_f,'b',t,kappa_f,'r')

 % Call of the function SEIQRDP.m with the fitted parameters
[S,E,I,Q,R,D,P] = SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,E0,I0,Q0,R0,D0,t,lambdaFun,kappaFun);
%------------------------------------------
%
% Comparison of fitted and real data
%
%------------------------------------------
figure
semilogy(time1,Q,'b',time1,R,'g');
hold on
semilogy(time(1:end-Nest),Active(1:end-Nest),'bo',time(1:end-Nest),Recovered(1:end-Nest),'go');
semilogy(time(length(Active)-Nest+1:length(Active)),Active(length(Active)-Nest+1:length(Active)),'b<',...
    time(length(Active)-Nest+1:length(Active)),Recovered(length(Active)-Nest+1:length(Active)),'g<');
% ylim([0,1.1*Npop])
ylabel('Number of cases')
xlabel('time (days)')
leg = {'Active (fitted)','Recovered (fitted)','Active (reported)','Recovered (reported)'};
legend(leg{:},'location','southoutside')
set(gcf,'color','w')
grid on
axis tight
set(gca,'yscale','lin')

figure
semilogy(time1,D,'b');
hold on
semilogy(time(1:end-Nest),Deaths(1:end-Nest),'bo');
semilogy(time(length(Active)-Nest+1:length(Active)),Deaths(length(Active)-Nest+1:length(Active)),'b<');
% ylim([0,1.1*Npop])
ylabel('Number of cases')
xlabel('time (days)')
leg = {'Deceased (fitted)','Deceased  (reported)'};
legend(leg{:},'location','southoutside')
set(gcf,'color','w')
grid on
axis tight
set(gca,'yscale','lin')
%tt=datenum(time(length(Active)-Nest+1:length(Active)))-datenum(time(length(Active)-Nest));
tt=datenum(time(1:length(Active)-Nest)-datenum(time(1)));

[C,i_t,i_tt]=intersect(t,tt);
re_R=sqrt(sum(R(i_t)-Recovered(i_tt)).^2./sum(Recovered(i_tt)).^2);
re=sqrt(sum(Q(i_t)-Active(i_tt)).^2./sum(Active(i_tt)).^2);
re_D=sqrt(sum(D(i_t)-Deaths(i_tt)).^2./sum(Deaths(i_tt)).^2);
RE=[RE;re];
RE_R=[RE_R;re_R];
RE_D=[RE_D;re_D];
checkRates(time(1:length(Active)-Nest),Active(1:length(Active)-Nest),Recovered(1:length(Active)-Nest),Deaths(1:length(Active)-Nest),kappaFun,lambdaFun,Kappa1,Lambda1);
end

%--------------------------------
%

figure (52)
semilogy(time1,E,'b',time1,I,'g');

ylabel('Number of cases')
xlabel('time (days)')
leg = {'Exposed (fitted)','Infectious (fitted)'};
legend(leg{:},'location','southoutside')
set(gcf,'color','w')
grid on
axis tight
set(gca,'yscale','lin')






