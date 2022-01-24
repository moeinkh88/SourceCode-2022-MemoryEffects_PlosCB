%% Figure S10
% BU BT
clear 
clc
%% Inputs
global A mu

order1=1:-.02:.84;
order2=1:-.02:.84;
X0= [.3; .37];
mu=[0.599 0.626];

t0=0;
T=10000;
h=.1;
F=@fun;
JF=@Jfun;

A=[-0.9059 -0.9377;-0.972 -0.9597];

%% fix points
xx1=[(A(1,2)*mu(2)-A(2,2)*mu(1))/(A(1,1)*A(2,2)-A(1,2)*A(2,1)),...
    (A(2,1)*mu(1)-A(1,1)*mu(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1))];
xx2=flip(xx1);
xx3=[0,-mu(2)/A(2,2)];
xx4=[-mu(1)/A(1,1),0];

x12=[xx1;xx2;xx3;xx4];

% Jac=JF(1,[x1,x2]);
% eig(Jac);

M1=length(order1);
M2=length(order2);
ConvergT=zeros(M1,M2);
for i=1:M1
    for j=1:M2
[t,X]=FDE_PI2_IM([order1(i),order2(j)],F,JF,t0,T,X0,h);

Err=braycd(X(:,end),x12');
[~,indFix]=min(Err);

indx=find(braycd(X,x12(indFix,:)')<5e-3);
ConvergT(i,j)=t(indx(1));
    end
end

%% plotting

figure
h=heatmap(1-order1,1-order2,ConvergT');
h.XLabel = 'Memory of BU';
h.YLabel = 'Memory of BT';

ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
