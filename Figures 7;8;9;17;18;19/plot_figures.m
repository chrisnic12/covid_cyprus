%%% Author Sergios Agapiou

%%% Code for plotting the results corresponding to Compartmental Model 1 in Section 2.4.1
%%% of "Modeling of Covid-19 Pandemic in Cyprus" by Agapiou et al. Results
%%% and plots contained in Section 3.5 and A-6


load('test.mat')


%t1=datetime(2020,3,4, 'Format', 'dd-MMM');
%t2=datetime(2020,5,26, 'Format', 'dd-MMM');
%%% For data on local transmission only need to change to
t1=datetime(2020,3,7, 'Format', 'dd-MMM');
t2=datetime(2020,5,29, 'Format', 'dd-MMM');
t=t1:t2;

close all

%%% Produces Figure 8 and 18
figure(1)
map = brewermap(num_periods,'Set1'); 
r=zeros(2,num_periods);
for j=1:num_periods
    r(:,j)=quantile(Re_vec(j,iter/5:end)',[0.025,0.975]);
    subplot(3,2,j)
    histogram(Re_vec(j,:),0:.2:7,'facecolor',map(j,:),'Normalization','probability')
    ylabel('Probability')
    axis tight
    ylim([0 0.31]);
    d1=(j-1)*num_times+1;
    d2=j*num_times;
    title(datestr(t(d1)) + " -- " + datestr(t(d2)))
end


d=datestr(t(1),'dd/mm') + " - " + datestr(t(num_times),'dd/mm');
for l=2:num_periods
    d=[d;datestr(t((l-1)*num_times+1),'dd/mm')+ " - " + datestr(t(l*num_times),'dd/mm')];
end

%%% Produces Figures 9 and 19
figure(2)
subplot(2,1,1)
y=median(Re_vec');
x=[1:num_periods];
errorbar([1:num_periods],y,r(1,:)-median(Re_vec'),r(2,:)-median(Re_vec'),'*','LineWidth', 1.5)
for k=1:numel(x)
      text(x(k),y(k),['  ', num2str(round(y(k),2))])
end
hold on
plot(0:11,ones(1,12),'--','LineWidth',1.3,'Color', 'k')
hold off
title('Median and $95\%$ Credible Intervals', 'Interpreter', 'latex','FontSize',14)
set(gca, 'XTick',[1 : num_periods],'XTickLabel',d);
set(gca,'XLim',[0.5 num_periods+.5],'YLim',[-0.2 max(r(2,:))+0.2])
xtickangle(30)
subplot(2,1,2)
prob=sum(Re_vec'<1)/length(Re_vec);
plot(prob,'+-', 'LineWidth', 1.7)
for k=1:numel(x)
      text(x(k)+0.05,prob(k)+0.05,[num2str(round(prob(k),2))])
end
set(gca,'XLim',[0.5 num_periods+.5],'YLim',[-0.05 1.05])
set(gca, 'XTick',[1 : num_periods],'XTickLabel',d);
xtickangle(30)
title('Posterior probability of $R_t<1$', 'Interpreter', 'latex','FontSize',14)


%%% Produces Figures 7 and 17
figure(3)
for j=1:num_periods
    subplot(2,3,j)

    x = 0:0.1:1; 
    y = 0:0.2:2; 
    [X,Y] = meshgrid(x,y); 

    pdf = hist3([alpha_vec(j,find(beta_vec(j,:)<2))', beta_vec(j,find(beta_vec(j,:)<2))'],{x y}); 
    pdf_normalize = (pdf'./length(alpha_vec(j,find(beta_vec(j,:)<2))')); 
    surf(X,Y,pdf_normalize) 
    d1=(j-1)*num_times+1;
    d2=j*num_times;
    caxis([0 0.085])
    title(datestr(t(d1)) + " -- " + datestr(t(d2)))
    xlabel('$\alpha$','Interpreter', 'latex')
    ylabel('$\beta$ (days$^{-1}$)','Interpreter', 'latex')
    colorbar
    view(2)
end