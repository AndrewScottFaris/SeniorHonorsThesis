% run_constrained_ramsey_simu returns 
%  policy functions and  simulated time series  under alternative equilibrium selection criteria and the Ramsey equilibrium in the economy presented in
% ``Multiple Equilibria in Open Economies with  Collateral Constraints,''  
% by Stephanie Schmitt-Grohé and Martín Uribe. 
%The output is in 12 .mat files. 
%©  Stephanie Schmitt-Grohé and Martín Uribe, January 2018.

clc
%help run_constrained_ramsey_simu
%disp('  ')
%disp('This code took about 4.5 hours to complete in 2019.'); disp(' ')
 
%Unpercentage the following line:
!copy tpm13.mat tpm.mat
%In this case the code takes 7 minutes to complete. 



%Economy with xi=0.5
disp('The Economy with xi=0.5')

disp('Computing Policy functions under equilibrium selection criteria, b, c, and ab (criterion ac is the same as ab in this calibration)')
Andrewconstrained_152018('b');disp('criterion b completed')
Andrewconstrained_152018('c');disp('criterion c completed')
Andrewconstrained_152018('ab');disp('criterion ab completed')

%disp('Computing Policy functions under the Ramsey equilibrium')
%ramsey_1122018; disp('ramsey equilibrium completed')

disp('Simulating time series')
simu_4192018('Andrewconstrained_152018b')
simu_4192018('Andrewconstrained_152018c')
simu_4192018('Andrewconstrained_152018ab') 
%simu_4192018('ramsey_1122018')

%Economy with xi=0.83
%disp('  ')
%disp('The Economy with xi=0.83')
%disp('Computing the policy function in the xi=0.83 economy under equilibrium selection criterium ab (ac, b, and c, is the same as ab in this calibration)')
%constrained_152018('ab',0.83,0.27,1.15); disp('criterion ab in xi=0.83 economy  completed')

%disp('Computing Policy functions under the Ramsey equilibrium')
%ramsey_1122018(0.83,0.27,1.15); disp('ramsey equilibrium in xi=0.83 economy completed')

%disp('Simulating time series when xi=0.83')
%simu_4192018('constrained_152018ab_xi83')
%simu_4192018('ramsey_1122018_xi83')


