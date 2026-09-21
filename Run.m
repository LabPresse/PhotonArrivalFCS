clc
clear
close all

Diffusion          =    10;
Lambda_P_real      =    100000  ;    % Emission rate per second
Lambda_b_real      =    1000    ;    % Background rate per second 
Rx                 =    0.5     ;    % Radious of space in micro
Ry                 =    Rx      ;
Rz                 =    2       ;
wxy                =    0.3    ;    % Semi-axis of confocal in micro
wz                 =    1.535   ;    % Semi-axis of confocal in micro in z direction
Number_particles   =    4      ;    % Number of particle in the space
Length_signal      =    5000    ;    % Number of steps
rc                 =    wxy     ;    % Semi-axis of confocal in micro
zc                 =    wz      ;    % Semi-axis of confocal in micro in z direction

minn               =  0.001     ;
maxx               =  0.05      ;
 
 
% [ Time , emittion_rate , X0, XX , Y0 , YY , Z0 , ZZ] = Sample_Generator( ...

[ Time , emission_rate , DeltaT , X0, XX , Y0 , YY , Z0 , ZZ] = Sample_Generator( ...
...
...
minn , maxx, Diffusion , Lambda_P_real , Lambda_b_real , Rx , Ry , Rz , wxy , ...
wz , Number_particles , Length_signal );
 


 

SSignal_time=Time(min(find(Time>=minn))+find(Time(find(Time>=minn))<=maxx)-1);
% DeltaT=(diff(SSignal_time))';

MINn=find(min(SSignal_time)==Time);
MAXn=find(max(SSignal_time)==Time);

XxX=XX(MINn+1:MAXn,:);
YyY=YY(MINn+1:MAXn,:);
ZzZ=ZZ(MINn+1:MAXn,:);

% Plot the Simulated data
  subplot(2,1,1) 
 for i=1:Number_particles
     plot(Time(2:end),emission_rate(2:end,i),'.-')
     hold on
 end
 xlabel('Time (s)', 'Fontsize',25)
 ylabel('Emission', 'Fontsize',25)
% set(gca,'XTick',[0 0.01 0.03 0.05])
    set(gca,'YTick',[0 1])
 set(gca,'FontSize',15)
 
subplot(2,1,2) 
stem(minn+[cumsum(DeltaT)],[ones(1,length(DeltaT))])
xlabel('Time (s)','Fontsize',15)
ylabel('Photon arrivals', 'Fontsize',15)
%     set(gca,'XTick',[0 0.01 0.03 0.05])
%     set(gca,'YTick',[0 1])
    set(gca,'FontSize',15)

%% Runing the Gibbs sampling
 
% It contains functions for FFBS , sampler for D, and x


number_species = 10                          ;


Prior_D_alpha       =    1                   ;               % Alpha for gamma prob. dist. for prior D 
Prior_D_beta        =    10                  ;               % Beta for gamma prob. dist. for prior  D
    
D                   =    1                   ;              % Initial Value for D


Lambda_P            =   Lambda_P_real        ;   % Initial Value for Lambda_P
Lambda_b            =   Lambda_b_real        ;   % Initial Value for Lambda_b
%   
b                   = ones(number_species,1) ;

q                   = ones(number_species,1) ;

alpha_q             =    1                   ;  
beta_q              =    1                   ;
%
M0                  =    1                   ;
var0                =    1                   ;
%
Prior_Lambda_P_alpha =   1                   ;
Prior_Lambda_P_beta  =   10000               ;
Lambda_P_alpha       =   10                  ;
% 
Prior_Lambda_b_alpha =   10                  ;
Prior_Lambda_b_beta  =   1000                ;
Lambda_b_alpha       =   200                 ;

x       = randn(number_species,length(DeltaT));
y       = randn(number_species,length(DeltaT));
z       = randn(number_species,length(DeltaT));

concen_radious    = [0.5 1 2 5] ;
concentration     = zeros(1,length(DeltaT),length(concen_radious)) ;


%%
tic

Max_iter             = 1000                 ;         % Number of iteration 

[D         , x                        , y                 , z                        , ...
 Lambda_P  , acceptance_rate_lambdaP  , Lambda_b          , acceptance_rate_lambda_b , ...
 b         , q                        , acceptance_rate_b , concentration          ] = ...
...
...
Gibbs_sampling_function_withNoise_Continous2( ...
...
...
DeltaT                ,  Max_iter               ,  Prior_D_alpha   ,  Prior_D_beta          ,  ...
D                     ,  Lambda_P               ,  Lambda_b        ,  rc                    ,  ...
zc                    ,  number_species         ,  alpha_q         ,  beta_q                ,  ...
b                     ,  q                      ,  M0              ,  var0                  ,  ... 
Prior_Lambda_P_alpha  ,  Prior_Lambda_P_beta    ,  Lambda_P_alpha  ,  Prior_Lambda_b_alpha  ,  ...
Prior_Lambda_b_beta   ,  Lambda_b_alpha         ,  x               ,  y                     ,  ...
z                     ,  concentration          ,  concen_radious        );

toc
%%
concen_radious    = [0.5 1 2] ;
% Measure the size of the data
sign_siz =  length(DeltaT) ;

% Create an empty matrix of concentration real
concentrationrel = zeros(length(concen_radious),sign_siz);

% Calculate the real concentration of particles 
for n=1:size(XX,2)
    for l=1:length(concen_radious)
        concentrationrel(l,:) = concentrationrel(l,:)+...
                                        ( sqrt( ((XxX(:,n)'./rc).^2)+...
                                                ((YyY(:,n)'./rc).^2)+...
                                                ((ZzZ(:,n)'./zc ).^2) )<=concen_radious(l) )...
                                       ./((pi^(3/2))*wxy*wxy*wz*(concen_radious(l)^3)*(10^-18)*(6.022*10^23)*(10^-6));
    end
end


mean25   =  [];
mean50   =  [];
mean75   =  [];

test=nan(size(concentration,1),1);

for l = 1 : length(concen_radious)
    for k = 1 : sign_siz
        mean25(l,k)=quantile(concentration(:,k,l)/((pi^(3/2))*wxy*wxy*wz*(concen_radious(l)^3)*(10^-18)*(6.022*10^23)*(10^-6)),0.25)      ;
        mean50(l,k)=quantile(concentration(:,k,l)/((pi^(3/2))*wxy*wxy*wz*(concen_radious(l)^3)*(10^-18)*(6.022*10^23)*(10^-6)),0.5)     ;
        mean75(l,k)=quantile(concentration(:,k,l)/((pi^(3/2))*wxy*wxy*wz*(concen_radious(l)^3)*(10^-18)*(6.022*10^23)*(10^-6)),0.75)     ;
    end
end

% Import the ell unnicode to the matlab
ell = char(hex2dec(strsplit('2113')));

% plot the real a learned concentration based on nano mol(nM) in different
% standard effective sizes of the confocal
for l = 1 : length(concen_radious)
    subplot(length(concen_radious)+1,1,l)
    
    plot(cumsum(DeltaT),concentrationrel(l,:),'g');
    hold on
    plot(cumsum(DeltaT),mean25(l,:),'--','color','r')
    hold on
    plot(cumsum(DeltaT),mean50(l,:),'b')
    hold on
    plot(cumsum(DeltaT),mean75(l,:),'--','color','r')
  
    xlabel('Time step','Fontsize',15)
    ylabel({'concentration','nM',[ell,'=',num2str(concen_radious(l))]},'Fontsize',17)
    box off
    set(gca,'XTick',[0 0.01 0.03 0.05])
    legend('Exact value of concentration','25-75% of posterior','Median of posterior')
%     set(gca,'YTick',[])
    set(gca,'FontSize',18)
end


subplot(length(concen_radious)+1,1,length(concen_radious)+1)
edge=cumsum([0,repmat(SSignal_time(end)/1000,1,1000)]);
[signalc,edges] = histcounts(SSignal_time,edge);
plot(edge(2:end),signalc)
% xlim([0 sign_siz])
% ylim([0 max(signal)+1])
% xlabel(['Time step (',num2str((SSignal_time(end)/1000)*1000000),'\mus)'],'Fontsize',15)
xlabel('Time step','Fontsize',15)
ylabel('Observed photons','Fontsize',15)
set(gca,'XTick',[0 0.01 0.03 0.05])
set(gca,'FontSize',18)




%%
subplot(2,1,1)
plot(sum(q))

subplot(2,1,2)
histogram(sum(q))

%% Load
subplot(2,1,1)
plot(sum(b))
xlabel('No. of Iterations')

subplot(2,1,2)
histogram(sum(b))
xlabel('No. of Particles')

%% Display the graphs
autocorr(D,100000)

%% D
% subplot(1,2,1)
plot(D);
%%
% hold on
% 
% plot([0 Max_iter] ,[Diffusion Diffusion])
% % axis([0 Max_iter 0 10])
% xlabel('No. of Iterations')
% ylabel('D [\mum^2/s]')
% 
% subplot(1,2,2)
% edges = [min(D):max(D)];
% histogram(D,40000, 'Normalization','Probability');
d_bnd = logspace(-1,2,100);
histogram(D,d_bnd,'Normalization','pdf')
set(gca,'Xscale','log')
% axis([0 500 0 .00065])
% set(gca,'Xscale')
hold on
% DD = 0:0.01:20;
% PriorD = 0.05.*(DD.^(-Prior_D_alpha-1)).*(exp(-Prior_D_beta./DD)).*(Prior_D_beta.^Prior_D_alpha)/gamma(Prior_D_alpha); 
% plot(DD,PriorD);
% hold on

y1=get(gca,'ylim');
plot([Diffusion Diffusion],y1,'--','LineWidth',3.5,'color','r')

xlabel('Diffusion coefficient (\mum^2/s)','Fontsize',22)
ylabel('Posterior Prob. Dist.','Fontsize',22)
% set(gca,'YTick',[0 0.005 0.01 0.015 0.02 0.025])
% set(gca,'YTick',[0 0.2 0.4 0.6])
% set(gca,'XTick',[1 10 100])
set(gca,'FontSize',18);
box off
hk =legend('Posterior using D=10 \mum^2/s','Exact Value of Diffusion Coef.');
legend('boxoff')

%% LAMBDA_P
% subplot(1,2,1)
% plot(Lambda_P); 
% % hold on
% % TTT=1:Max_iter; Lambda_P_real1(TTT)=Lambda_P_real; plot(Lambda_P_real1)
% xlabel('No. of Iterations')
% ylabel('\lambda_P [#/s]')
% 
% subplot(1,2,2)
% histogram(5*Lambda_P(:,[floor(0.3*length(Lambda_P)):end]),500,'Normalization','Probability'); 
d_bnd = logspace(3,7,100);
histogram(Lambda_P(:,[floor(0.3*length(Lambda_P)):end]),d_bnd,'Normalization','pdf')
hold on
xlabel('\lambda_P')
set(gca,'Xscale','log')
% hold on
% DD = 0:0.01:20;
% PriorD = 0.05.*(DD.^(-Prior_D_alpha-1)).*(exp(-Prior_D_beta./DD)).*(Prior_D_beta.^Prior_D_alpha)/gamma(Prior_D_alpha); 
% plot(DD,PriorD);
hold on

y2=get(gca,'ylim');
plot([Lambda_P_real Lambda_P_real],y2,'--','LineWidth',3.5,'color','r')
xlabel('Photon Emission Rate (photon/s)','Fontsize',15)
ylabel('Posterior Prob. Dist.','Fontsize',15)
set(gca,'FontSize',18);
set(gca,'XTick',[10^4 10^5 10^6 10^7])
% set(gca,'YTick',10^-5*[0.1 0.2 0.3 0.4 0.5 0.6])
box off
hk1 =legend('Posterior using \mu_{part}=10^5', 'Exact value of Emission Rate');
legend('boxoff')

%% LAMBDA_b
subplot(1,2,1)
plot(Lambda_b); hold on
% TT=1:Max_iter; Lambda_b_real1(TT)=Lambda_b_real; plot(Lambda_b_real1)
xlabel('No. of Iterations')
ylabel('\lambda_b [1/s]')

subplot(1,2,2)
histogram(Lambda_b(:,[0.3*length(Lambda_b):end]),100,'Normalization','Probability'); hold on
xlabel('\lambda_b')
hold on

y3=get(gca,'ylim');
plot([Lambda_b_real Lambda_b_real],y3,'--','LineWidth',1.75,'color','r')
xlabel('Bachground rate (1/s)','Fontsize',15)
ylabel('Posterior','Fontsize',15)
box off
hk2 =legend('Posterior prob. dist.');

