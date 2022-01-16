clear 
clc

global n N Ki b Kij
%% Coefficients and Conditions

N=2;

p=1:-.1:.6;

n=2; % Hill coefficient

Kij=0.1*ones(2); % interaction matrix

Ki=1*ones(N,1); % death rate

T=300; %  final time

b=[1, 2]; % growth rates for cases: False, Pulse, and Periodic

% X0=[.9;.15]; % initial conditions
X0=[.8;.2];

t0=0;
h=.2;
F=@funGonze;

for j=1:length(p)
    order=p(j)*ones(2,1);
[t,x]=FDE_PI12_PC(order,F,t0,T,X0,h);
X(j,:,:)=x;
end

%%
figure

x1(:,:)=X(1,:,:)./sum(X(1,:,:));
x2(:,:)=X(2,:,:)./sum(X(2,:,:));
x3(:,:)=X(3,:,:)./sum(X(3,:,:));
x4(:,:)=X(4,:,:)./sum(X(4,:,:));
x5(:,:)=X(5,:,:)./sum(X(5,:,:));
hold on
p1=plot(t,x1(1,:),'Color',[0,0,.7],'LineWidth',3,'DisplayName','0');
p2=plot(t,x2(1,:),'Color',[0,0,.9],'LineWidth',3,'DisplayName','0.1');
p3=plot(t,x3(1,:),'Color',[0,0.5,1],'LineWidth',3,'DisplayName','0.2');
p4=plot(t,x4(1,:),'Color',[0,.7,1],'LineWidth',3,'DisplayName','0.3');
p5=plot(t,x5(1,:),'Color',[0,.9,1],'LineWidth',3,'DisplayName','0.4');

p6=plot(t,x1(2,:),'Color',[.4,0,0],'LineWidth',3,'DisplayName','0');
p7=plot(t,x2(2,:),'Color',[.6,0,0],'LineWidth',3,'DisplayName','0.1');
p8=plot(t,x3(2,:),'Color',[.8,0,0],'LineWidth',3,'DisplayName','0.2');
p9=plot(t,x4(2,:),'Color',[1,0,0],'LineWidth',3,'DisplayName','0.3'); 
p10=plot(t,x5(2,:),'Color',[1,.3,0],'LineWidth',3,'DisplayName','0.4');

set(gca,'FontSize',14)



xlabel('Time')
ylabel('Relative abundance')

set(gca,'FontSize',14)
% legend('0','0.1','0.2', '0.3', '0.4')

leg = legend('show');
title(leg,'X1    \mu   X2   \mu ')
leg.NumColumns = 2;