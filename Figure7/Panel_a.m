clc
clear

F=@fun1;
al=1:-.1:.6;
h=0.01;

for i=1:length(al)
    [t,x] = FDE_PI12_PC(al(i),F,0,10,.1,h);
    X(i,:)=x;
end
%%
figure
p1=plot(t,X);
xlabel('Time')
ylabel('Logistic growth curve')
set(p1,'LineWidth',3)
set(gca,'FontSize',14)
legend('0','0.1','0.2', '0.3', '0.4')

leg = legend('show');
title(leg,'Memory')
