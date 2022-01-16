%% Main_ThreeSpecies
%        ------------------------------------------------------------------
%                   This code solves a there species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=X_i(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                   
%        ------------------------------------------------------------------
%        mu - Order of derivatives, [mu_B,mu_R,mu_G]  0<mu(i)=<1, e.g. mu=[1,.2,1];
%        ------------------------------------------------------------------
%        n -  Hill coefficient, e.g. n=2;
%        ------------------------------------------------------------------
%        N -  Number of Species, e.g. N=3;
%        ------------------------------------------------------------------
%        Kij - Interation matrix, e.g. Kij=0.1*ones(N);
%        ------------------------------------------------------------------
%        Ki - Death rate, e.g. Ki=1*ones(N,1);
%        ------------------------------------------------------------------
%        T - Final time, e.g. T=600;
%        ------------------------------------------------------------------
%        x0 - Initial conditions, e.g. x0=[1/3;1/3;1/3];
%        ------------------------------------------------------------------
%        Perturbation - This is changes of the species growth rates, e.g. Perturbation='OUP';
%                       Possible usages 'False', 'Pulse, 'Periodic', 'OUP', and 'OUP_new'                      
%                       False: No perturbation
%                       Pulse1: A pulse in (20,60) for Figure 2a
%                       Pulse2: Similar to Pulse1 with a greater strength for Figure 2b
%                       Pulse3: Two pulses in (60,100) and (200,330) for Figure 3a
%                       Pulse4: Two pulses in (60,100) and (400,530) for Figure S2a
%                       Periodic: Periodic perturbation with 20 span for Figure 3b
%                       OUP: Stochastic pertubation used in the paper; requirement: T=<700, For Figure 4b-c & S2b                          
%                       OUP_new: New generating stochastic perturbation
%        ------------------------------------------------------------------
%        b - Growth rates for cases: False, Pulse, and Periodic, e.g. b=[1, .95, 1.05];
%
%---------------------------------------
% Outputs
%        t - Simulated time interval
%        x - Species abundances 
%        B - Growth rates including perturbation
%
%
%  Please, report any problem or comment to :
%          moein dot khalighi at utu dot fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc
global n N Ki b Kij T x0
%% Coefficients and Conditions

mu=0.9*[1,1,1]; % [mu_B,mu_R,mu_G] Order of derivatives,  0<mu(i)=<1

n=2; % Hill coefficient

N=3; % number of species

Kij=0.1*ones(N); % interaction matrix

Ki=1*ones(N,1); % death rate

T=700; %  final time

Perturbation='OUP'; % Possible usages 'False', 'Pulse1', 'Pulse2', 'Pulse3', 'Pulse4', 'Periodic', 'OUP', and 'OUP_new'

b=[1, .95, 1.05]; % growth rates for cases: False, Pulse, and Periodic

x0=1/3*[1,1,1]'; % initial conditions

Fun=@fun3;
        clear b        
        load('b3OUP.mat');
        bb=@(t,N)b(t,N);
        B=b;
        
        t0=0; % initial time
h=0.01; % step size for computing


% solver for fractional differential equation
[t, x] = FDE_PI12_PC(mu,Fun,t0,T,x0,h);

%% plotting

RelX=x./(ones(N,1)*sum(x)); % Relative abundances

f=figure;
f.Renderer='painters';

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
        
    %%plot of growth rate with stochastic perturbation
    
    tt=1:T; % time simulating with 1 step size  
    
    % Plot the simulated growth rates
    
    pb=plot(tt,B(1:T,1),'b',tt,B(1:T,2),'r',tt,B(1:T,3),'g');
    
    % Settings for plot
    set(pb,'LineWidth',4)
    ylabel('Growth Rates: OUP')
    xlabel('Time')
    legend('b_B','b_R', 'b_G')
    set(gca,'FontSize',23, 'FontWeight', 'bold')
    
    %%plotting relative abundance of species
       figure
            p=plot(t,RelX(1,:),'b',t,RelX(2,:),'r',t,RelX(3,:),'g');
       ylabel('Relative abundance')


% Settings for the general plot
set(p,'LineWidth',5)
set(gca,'FontSize',23, 'FontWeight', 'bold')
xlabel('Time')

% % Generate legend for showing memory of each species
legend(['X_{B}/X_{total}, Memory_B=',num2str(1-mu(1))],...
    ['X_{R}/X_{total}, Memory_R=',num2str(1-mu(2))],...
    ['X_{G}/X_{total}, Memory_G=',num2str(1-mu(3))],'Location', 'Best');



