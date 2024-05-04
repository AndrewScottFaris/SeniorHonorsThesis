%self_fulfilling_b.m produces Figure 9  (Typical Self-Fulfilling Crises: Equilibrium Selection Criterion (C)) in the paper  
%``Multiple Equilibria in Open Economies with  Collateral Constraints,''  
% by Stephanie Schmitt-Grohé and Martín Uribe. 
%©  Stephanie Schmitt-Grohé and Martín Uribe, January 2018.

%find self-fulfilling crises, that is, crises in C that are not crises inA

clear all
clf

rows = 4;cols=2; 

orient tall

load  simu_Andrewconstrained_152018c.mat CRISIS YT P CT CA Dp YN TB SLACK CA rstar D S CRISIS
%produced by running
%simu_4192018('constrained_152018c')
%In turn, the file constrained_152018c.mat is produced by running constrained_152018.m and setting eqm_selection_criterion = 'c';
v3=[-0.2 0.4];

ww = 4; %window

s1 = find(CRISIS==1);

find(s1<=ww|s1>=length(YT)-ww);
s1 = setdiff(s1,s1(ans));
nw = numel(s1);

WP1 = zeros(nw,2*ww+1);
WYT1 = zeros(nw,2*ww+1);
WYN1 = zeros(nw,2*ww+1);
WCT1 = zeros(nw,2*ww+1);
WCA1 = zeros(nw,2*ww+1);
WDp1 = zeros(nw,2*ww+1);
WTB1 = zeros(nw,2*ww+1);
WSLACK1 = zeros(nw,2*ww+1);
WCRISIS1 = zeros(nw,2*ww+1);
WD1 = zeros(nw,2*ww+1);
WS1 = zeros(nw,2*ww+1); %states

for i=1:nw
WP1(i,:) = P(s1(i)-ww:s1(i)+ww)';
WYT1(i,:) = YT(s1(i)-ww:s1(i)+ww)';
WYN1(i,:) = YN(s1(i)-ww:s1(i)+ww)';
WCT1(i,:) = CT(s1(i)-ww:s1(i)+ww)';
WCA1(i,:) = CA(s1(i)-ww:s1(i)+ww)';
WDp1(i,:) = Dp(s1(i)-ww:s1(i)+ww)';
WTB1(i,:) = TB(s1(i)-ww:s1(i)+ww)';
WSLACK1(i,:) = SLACK(s1(i)-ww:s1(i)+ww)';
WCRISIS1(i,:) = CRISIS(s1(i)-ww:s1(i)+ww)';
WD1(i,:) = D(s1(i)-ww:s1(i)+ww)';
WS1(i,:) = S(s1(i)-ww:s1(i)+ww)';
end


%dynamics under equilibrium selection criterion  A

AP1 = zeros(nw,2*ww+1);
AYT1 = zeros(nw,2*ww+1);
AYN1 = zeros(nw,2*ww+1);
ACT1 = zeros(nw,2*ww+1);
ACA1 = zeros(nw,2*ww+1);
ADp1 = zeros(nw,2*ww+1);
ATB1 = zeros(nw,2*ww+1);
ASLACK1 = zeros(nw,2*ww+1);
AD1 = zeros(nw,2*ww+1);
AS1 = zeros(nw,2*ww+1); %states
ACRISIS1 = zeros(nw,2*ww+1);

%Construct  time series for A for those windows 

clear d dp pXgridlevel yNgridlevel crisis p cT dgrid dpix ny slack ca tb slack 

load Andrewconstrained_152018ab.mat d dp pXgridlevel yNgridlevel crisis p cT ny  dgrid dpix slack rstar 
%produced by running
% constrained_152018.m and setting eqm_selection_criterion = 'ab' (for this calibration criterion ac is identical to ab)

%current account: 
%Construct the current account from the formula ca_t = 1/(1+r_t-1)d_t - 1/(1+r_t)d_t+1
lagg(Dp./(1+rstar));
[0;ans(:,2)];
CA = ans - Dp./(1+rstar);

for w=1:nw
yT0 = WYT1(w,1); 
yN0 = WYN1(w,1); 
d0 = WD1(w,1); 

    for t=1:2*ww+1; 
AYT1(w,t) = yT0;
AYN1(w,t) = yN0;
AD1(w,t) = d0;
i = find(pXgridlevel==yT0 & yNgridlevel==yN0);
tolerance = 1e-4;
j = find(abs(dgrid-d0)<tolerance);
ACRISIS1(w,t) = crisis(i,j);
AP1(w,t) = p(i,j);
ACT1(w,t) = cT(i,j);
ATB1(w,t) = yT0 - cT(i,j);
ASLACK1(w,t) = slack(i,j);
AS1(w,t) = (j-1)*ny+i;
ADp1(w,t) = dgrid(dpix(i,j));
ACA1(w,t) = AD1(w,t)/(1+rstar) - ADp1(w,t)/(1+rstar);
d0 = dp(i,j); 
if t<2*ww+1
yT0 = WYT1(w,t+1);
yN0 = WYN1(w,t+1); 
end %if t<2*ww+1
end %for t
end %for w=1:nw

%now find those windows in which under policy A there was not no crisis 
Aix = find(sum(ACRISIS1, 2)==0);

thick = 1;
lstyle = '-'; 

Wcollateral = WSLACK1+WDp1;

subplot(rows,cols,1)
x =  mean(WYT1(Aix,:));
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Price of Exportables, $p^X_t$');
set(t, 'Interpreter', 'Latex')
xlim([-ww ww])
set(gca, 'XTick',[-5:5])
v = [0.85 1.05];
ylim(v)
hold on

x =  mean(WYN1(Aix,:));
subplot(rows,cols,2)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Nontraded Output, $y^N_t$');
set(t, 'Interpreter', 'Latex'), xlim([-ww ww]) ,set(gca, 'XTick',[-5:5]) 
ylim(v)
hold on

x =  mean(WCT1(Aix,:));
subplot(rows,cols,3)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Consumption of Tradables, $c^T_t$');
set(t, 'Interpreter', 'Latex'), xlim([-ww ww]) ,set(gca, 'XTick',[-5:5]) 
hold on

x =  mean(WP1(Aix,:));
subplot(rows,cols,4)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Relative Price of Nontradables, $p^N_t$');
set(t, 'Interpreter', 'Latex'), xlim([-ww ww]) ,set(gca, 'XTick',[-5:5]) 
hold on

x =  mean(Wcollateral(Aix,:));
subplot(rows,cols,5)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Collateral, $\kappa (y^M_t + p^X_ty^X_t + p^N_t y^N_t)$');
set(t, 'Interpreter', 'Latex'), xlim([-ww ww])   ,set(gca, 'XTick',[-5:5]) 
hold on

x =  mean(WDp1(Aix,:)/(1+rstar));
subplot(rows,cols,6)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Debt, $d_{t+1}/(1+r)$');
set(t, 'Interpreter', 'Latex'), xlim([-ww ww]) ,set(gca, 'XTick',[-5:5]) 

hold on

x =  mean(WTB1(Aix,:));
subplot(rows,cols,7)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Trade Balance, $y^T_t-c^T_t$');
set(t, 'Interpreter', 'Latex'), xlim([-ww ww]) ,set(gca, 'XTick',[-5:5]) 
ylim(v3); 
hold on

x =  mean(WCA1(Aix,:));
subplot(rows,cols,8)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
t=title('Current Account, $(d_t- d_{t+1})/(1+r)$');
set(t, 'Interpreter', 'Latex'), xlim([-ww ww]) ,set(gca, 'XTick',[-5:5]) 
ylim(v3)
hold on

lstyle = '--'; 
Acollateral = ASLACK1+ADp1;

subplot(rows,cols,1)
x =  mean(AYT1(Aix,:));
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
ylim(v)

x =  mean(AYN1(Aix,:));
subplot(rows,cols,2)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);
ylim(v)

x =  mean(ACT1(Aix,:));
subplot(rows,cols,3)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);

x =  mean(AP1(Aix,:));
subplot(rows,cols,4)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);

x =  mean(Acollateral(Aix,:));
subplot(rows,cols,5)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);

x =  mean(ADp1(Aix,:)/(1+rstar));
subplot(rows,cols,6)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);

x =  mean(ATB1(Aix,:));
subplot(rows,cols,7)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);

x =  mean(ACA1(Aix,:));
subplot(rows,cols,8)
plot(-ww:ww,x,'linewidth',thick, 'linestyle', lstyle);


shg




