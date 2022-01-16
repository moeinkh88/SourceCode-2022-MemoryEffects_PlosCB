clear 
clc

global A mu

order1=1:-.02:.9;
order2=1:-.02:.9;

mu=[0.599 0.626];

t0=0;
T=300;
h=.1;
F=@fun;
JF=@Jfun;

A=[-0.9059 -0.9377;-0.972 -0.9597];

%% fix points
xx1=[(A(1,2)*mu(2)-A(2,2)*mu(1))/(A(1,1)*A(2,2)-A(1,2)*A(2,1));...
    (A(2,1)*mu(1)-A(1,1)*mu(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1))];
xx3=[1e-3;-mu(2)/A(2,2)];
X0=xx3; % initial conditions


M1=length(order1);
M2=length(order2);
Resistance=zeros(M1,M2);

for i=1:M1
    for j=1:M2
        tic
        p=.24; %perturb
while 1
[t,X]=FDE_PI2_IM([order1(i),order2(j)],F,JF,t0,T,X0,h,p);

Df=diff(X(:,end-1:end)');
% plot(t,X)
if Df(1)<=0 && Df(2)>=0
    p=p+0.001;
else
    Resistance(i,j)=p-0.001;
    break
end

end
toc
    end
end

%%
figure
h=heatmap(1-order1,1-order2,Resistance');
h.XLabel = 'Memory of BU';
h.YLabel = 'Memory of BT';

ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
