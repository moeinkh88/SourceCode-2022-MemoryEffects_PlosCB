clear
clc
global n N b Ki

load('bb.mat');b=bb;
load('x0.mat');
load('Kij.mat');
n=2; % Hill coefficient
N=15; % number of species

Blue=1:5;Red=6:10;Green=11:15;

Ki=1*ones(N,1); % death rate
t0=0; T=700; % t0 is start time, and T is the end time

alpha=ones(N,1); % order of derivatives
alpha(Blue)=1-.2;
 
h=0.01; % step size for computing

[tt, xx] = FDE_PI12_PC(alpha,@fun,t0,T,x0,h);


%%
PcG= [0.18,0.40,0.14];
PcR= [0.92,0.27,0.18];

t=tt(1:50:end);
x=xx(:,1:50:end);

f=figure;
f.Renderer='painters';
pb=semilogx(t,x(Blue,:),'b');
set(pb,'LineWidth',2.3)
hold on
pr=semilogx(t,x(Red,:),'Color',PcR);
set(pr,'LineWidth',2.3)
pg=semilogx(t,x(Green,:),'Color',PcG);
set(pg,'LineWidth',2.3)
line([t0 T],[0.5,0.5],'LineStyle','--', 'color', 'k')
hold off

axis([t0 T 0 1])
xlabel('{time}','FontSize',15)

ylabel('Abundnace','FontSize',15)

set(gca,'Fontsize',24)

for i=1:5
    if ismember(i,2:5)==1
set(get(get(pb(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(pr(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(pg(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
% legend({'X_B', 'X_R', 'X_G'},'FontSize',14)


xticks([10^-2 10^0 10^2])
xticklabels({'-2','0','2'})
         
Pos = [795 359 700 457];
set(0, 'DefaultFigurePosition', Pos);
text(.45,.7,['Memory_{B}=',num2str(1-alpha(1))],'FontSize',20)
text(.45,.6,['Memory_{R}=',num2str(1-alpha(6))],'FontSize',20)
text(.45,.5,['Memory_{G}=',num2str(1-alpha(15))],'FontSize',20)
