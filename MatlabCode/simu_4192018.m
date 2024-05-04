function  simu_4192018(filename,T,Tburn)
%simu_4192018.m  
%Produces  artificial time series of length T after burning Tburn periods, implied by equilibrium policy functions in filename. 
%This program was written for the economy studied in
%``Multiple Equilibria in Open Economies  with Collateral Constraints,''  
% by Stephanie Schmitt-Grohé and Martín Uribe. 
%©  Stephanie Schmitt-Grohé and Martín Uribe, January 2018.

%read policy functions
eval(['load ' filename]);
%policy function, filename.mat,  is produced by running constrained_152018.m   or by
%running  ramsey_1122018.m   


if ~exist('slack','var')
slack = zeros(ny,nd);
end

if ~exist('crisis','var')
crisis = zeros(ny,nd);
end

if ~exist('mu','var')
    mu = zeros(ny,nd);
end

rng(1); %set seed of random number generator 

if nargin<2
T = 1e6; %length of simulated time series
end 

if nargin<3
Tburn = round(T/10); %burn-in period
end

%Set initial conditions of tradable output,  interest rate, and debt
aux1=unique(pXgridlevel);aux2=numel(aux1);
yT0 = aux1(floor(aux2/2)); %take the middle value of the marginal grid of the level of yT (position 3 in a vector of 4 elements)
aux1=unique(yNgridlevel);aux2=numel(aux1);
yN0 = aux1(floor(aux2/2)); %take the middle value of the marginal grid  of the level of yN (position 9 in a vector of 16 elements)
d0 = dgrid(floor(nd/2));

%pai is the transition probability matrix of (yT,r). Cpai is the cumulative probability matrix (useful for drawing of exogenous states yT and r)
Cpai = cumsum(pai,2); 


%Initializations
S = zeros(T+Tburn,1); %state
LA = zeros(T+Tburn,1);%marginal utility of cT
YT = zeros(T+Tburn,1); %yT
D = zeros(T+Tburn,1); %current debt
Dp = zeros(T+Tburn,1);%next-period debt
CRISIS = zeros(T+Tburn,1); %1 if near binding
P = zeros(T+Tburn,1); %relative price nontradables 
CT = zeros(T+Tburn,1);  %consumption of tradables 
C = zeros(T+Tburn,1); %composite consumption
YN = zeros(T+Tburn,1); %Country interest rate

SLACK = zeros(T+Tburn,1); %slackness of collateral constraint
MU = zeros(T+Tburn,1); %capital control tax rate

%Construct simulated time series
for t=1:T+Tburn
YT(t,1) = yT0;
YN(t,1) = yN0;
D(t,1) = d0;
i = find(pXgridlevel==yT0 & yNgridlevel==yN0);
tolerance = 1e-5;
j = find(abs(dgrid - d0) < tolerance);
CRISIS(t,1) = crisis(i,j);
P(t,1) = p(i,j);
CT(t,1) = cT(i,j);
C(t,1) = c(i,j);
LA(t,1) = la(i,j);

MU(t,1) = mu(i,j);
SLACK(t,1) = slack(i,j);
S(t,1) = (j-1)*ny+i;
Dp(t,1) = dgrid(dpix(i,j));
d0 = dp(i,j); 
%draw realization of exogenous state conditional on current exogenous state
L = sum(Cpai(i,:)<rand)+1;
yT0 = pXgridlevel(L);
yN0 = yNgridlevel(L);
end %for t

initial_state = S(1);

%Construct the current account from the formula ca_t = 1/(1+r_t-1)d_t - 1/(1+r_t)d_t+1
lagg(Dp./(1+rstar));
[0;ans(:,2)];
CA = ans - Dp./(1+rstar);

%Eliminate burn-in periods
YT = YT(Tburn+1:end);
D = D(Tburn+1:end);
Dp = Dp(Tburn+1:end);

P = P(Tburn+1:end);
CRISIS = CRISIS(Tburn+1:end);
C = C(Tburn+1:end);
CT = CT(Tburn+1:end);
LA = LA(Tburn+1:end);

SLACK = SLACK(Tburn+1:end);
S = S(Tburn+1:end);
YN = YN(Tburn+1:end);
CA = CA(Tburn+1:end);
MU = MU(Tburn+1:end);

%Produce simulations of other   variables of the model
TB = YT - CT; %trade balance
Y = YT + YN.*P; %output measured in terms of tradables 
DY = D./Y; %debt-to-output ratio (annual)
EDY = mean(DY); %mean debt-to-output ratio
TBY = TB ./Y; %trade balance to output ratio
 

%Unconditional probability distribution of debt
for i=1:nd
lad(i,1) = mean(D==dgrid(i)); %prob. distribution of debt
end

eval(['save simu_' filename '.mat']);

