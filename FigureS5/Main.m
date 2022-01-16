clear
clc
global n N b Ki Kij 

n=2; % Hill coefficient
N=15; % number of species

Kij=.2*ones(N); 
% Blue=1:10;Red=11:20;Green=21:30;
Blue=1:5;Red=6:10;Green=11:15;

%%%
Ki=.5*ones(N,1); % death rate
t0=0; T=1000; % t0 is start time, and T is the end time

s = RandStream('twister','Seed',123);
    RandStream.setGlobalStream(s)
    
b=.1.*abs(randn(1,N)) + 1; % growth rate

x0=ones(N,1)/4; % initial conditions (1/4)

alpha=ones(N,1); % order of derivatives

h=0.05; % step size for computing

[t, x] = FDE_PI12_PC(alpha,@fun,t0,T,x0,h);


%%

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];

for i=1:N
    Kij(i,i)=0;
end

xx=x./(ones(N,1)*sum(x));
% xx=x;

figure

subplot(1,3,1)
% pb=semilogx(t,x(Blue,:),'b');
pb=plot(t(1:20:end),xx(Blue,1:20:end),'b');
set(pb,'LineWidth',2)
hold on
% pr=semilogx(t,xx(Red,:));
pr=plot(t(1:20:end),xx(Red,1:20:end));
set(pr,'LineWidth',2,'color',PcR)
% pg=semilogx(t,xx(Green,:));
pg=plot(t(1:20:end),xx(Green,1:20:end));
set(pg,'LineWidth',2,'color',PcG)
% line([t0 T],[0.5,0.5],'LineStyle','--', 'color', 'k')
hold off
xlabel('Time','FontSize',15)
ylabel('Abundance','FontSize',15) 
% title(['Memory_{B}=',num2str(1-alpha(1)),', Memory_{R}=',num2str(1-alpha(6)),...
%     ', Memory_{G}=',num2str(1-alpha(15)),...
%     '\newline X_{B}(0)=',num2str(x0(1)),', X_{R}(0)=',num2str(x0(6)),...
%     ', X_{G}(0)=',num2str(x0(11))],'FontSize',14)
title(['Memory_{B}=',num2str(1-alpha(1)),', Memory_{R}=',num2str(1-alpha(6)),...
    ', Memory_{G}=',num2str(1-alpha(15))],'FontSize',14)
set(gca,'Fontsize',15)

set(gcf,'renderer','Painters')

axis tight
%%
alpha(Blue)=1; alpha(Red)=.4;alpha(Green)=.7;

h=0.05; % step size for computing

[t, x] = FDE_PI12_PC(alpha,@fun,t0,T,x0,h);


%%

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];

for i=1:N
    Kij(i,i)=0;
end

xx=x./(ones(N,1)*sum(x));
% xx=x;

subplot(1,3,2)
% pb=semilogx(t,x(Blue,:),'b');
pb=plot(t(1:20:end),xx(Blue,1:20:end),'b');
set(pb,'LineWidth',2)
hold on
% pr=semilogx(t,xx(Red,:));
pr=plot(t(1:20:end),xx(Red,1:20:end));
set(pr,'LineWidth',2,'color',PcR)
% pg=semilogx(t,xx(Green,:));
pg=plot(t(1:20:end),xx(Green,1:20:end));
set(pg,'LineWidth',2,'color',PcG)
% line([t0 T],[0.5,0.5],'LineStyle','--', 'color', 'k')
hold off
xlabel('Time','FontSize',15)
ylabel('Abundance','FontSize',15) 
% title(['Memory_{B}=',num2str(1-alpha(1)),', Memory_{R}=',num2str(1-alpha(6)),...
%     ', Memory_{G}=',num2str(1-alpha(15)),...
%     '\newline X_{B}(0)=',num2str(x0(1)),', X_{R}(0)=',num2str(x0(6)),...
%     ', X_{G}(0)=',num2str(x0(11))],'FontSize',14)
title(['Memory_{B}=',num2str(1-alpha(1)),', Memory_{R}=',num2str(1-alpha(6)),...
    ', Memory_{G}=',num2str(1-alpha(15))],'FontSize',14)
set(gca,'Fontsize',15)

set(gcf,'renderer','Painters')

axis tight
%%
alpha(Blue)=.5; alpha(Red)=.6;alpha(Green)=1;

h=0.05; % step size for computing

[t, x] = FDE_PI12_PC(alpha,@fun,t0,T,x0,h);


%%

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];

for i=1:N
    Kij(i,i)=0;
end

xx=x./(ones(N,1)*sum(x));
% xx=x;

subplot(1,3,3)
% pb=semilogx(t,x(Blue,:),'b');
pb=plot(t(1:20:end),xx(Blue,1:20:end),'b');
set(pb,'LineWidth',2)
hold on
% pr=semilogx(t,xx(Red,:));
pr=plot(t(1:20:end),xx(Red,1:20:end));
set(pr,'LineWidth',2,'color',PcR)
% pg=semilogx(t,xx(Green,:));
pg=plot(t(1:20:end),xx(Green,1:20:end));
set(pg,'LineWidth',2,'color',PcG)
% line([t0 T],[0.5,0.5],'LineStyle','--', 'color', 'k')
hold off
xlabel('Time','FontSize',15)
ylabel('Abundance','FontSize',15) 
% title(['Memory_{B}=',num2str(1-alpha(1)),', Memory_{R}=',num2str(1-alpha(6)),...
%     ', Memory_{G}=',num2str(1-alpha(15)),...
%     '\newline X_{B}(0)=',num2str(x0(1)),', X_{R}(0)=',num2str(x0(6)),...
%     ', X_{G}(0)=',num2str(x0(11))],'FontSize',14)
title(['Memory_{B}=',num2str(1-alpha(1)),', Memory_{R}=',num2str(1-alpha(6)),...
    ', Memory_{G}=',num2str(1-alpha(15))],'FontSize',14)
set(gca,'Fontsize',15)

set(gcf,'renderer','Painters')

axis tight