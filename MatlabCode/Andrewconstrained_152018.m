function Andrewconstrained_152018(eqm_selection_criterion,xi,dlower,dupper, xiMX, yXendowmentconstant)
%constrained_152018(eqm_selection_criterion,xi,dlower,dupper) 
%uses a policy-function iteration procedure to approximate the equilibrium of  the open economy with a flow 
%collateral constraint presented in 
%``Multiple Equilibria in Open Economies with  Collateral Constraints,''  
% by Stephanie Schmitt-Grohé and Martín Uribe. 
%Inputs
%eqm_selection_criterion to pick an equilibrium selection criterion, set the  string variable  follows:
%'b' favors the least severe of the two constrained outcomes (point B in the graphs shown in the paper cited above). 
%'c' favors the most severe of the two constrained outcomes (point C in the graphs) 
%'ab' and 'ac' favor unconstrained outcomes (point A in the graphs of the cited paper).  
% 'ab' is a criterion under which if at a given current state  the equilibrium conditions cannot be satisfied with a nonbinding collateral constraint but all equilibrium conditions are satisfied at two binding points, then criterion 'ab'  favor the least severe of the two constrained outcomes (point B in the graphs) . 
%Criterion 'ac' is defined in a similar fashion. Criteria 'ab' and 'ac' deliver the same equilibrium in the calibration considered here. 
%xi intratemporal elasticity of substitution of tradables for nontradables.
%dlower lower bound of debt grid.
%dupper upper bound of debt grid.
%tpmname transition probability matrix of exogenous state
%©  Stephanie Schmitt-Grohé and Martín Uribe, January 2018.

format compact

filename = ['Andrewconstrained_152018' eqm_selection_criterion];

load tpm.mat yTgrid yNgrid pai
pXgrid = yTgrid;

%%IDEA: we have a joint process of terms of trade: pXt and nontradable
ny = numel(pXgrid); %number of grid points for the joint process of tradable output and  nontradable output

%calibration:
rstar = 0.04; %interest rate
betta = 0.91; %subjective discount factor 
a = 0.05; %weight on traded consumption in the CES aggregator
b = 0.5; %weight on importable consumption in the CES aggregator
sigg = 2; %intertemporal elasticity of consumption substitution
kapa = 0.32*(1+rstar); %fraction of output that is pledgeable as collateral
yXendowmentconstant = 1; 
thet = 1; %importable good endowment



%Grids in levels
pXgridlevel = exp(pXgrid(:)); %grid of level of terms of trade
yNgridlevel = exp(yNgrid(:)); %grid of level of country interest rate

%natural debt limit
NdL = min(pXgridlevel)*(1+rstar)/rstar; 

%The debt grid: 
nd = 800;  %number of points

if nargin<2 | isempty(xi)
xi = 0.5;  %Elasticity of substitution between traded and nontraded goods
end


if nargin<3|isempty(dlower)
dlower =0.2;  %lower bound of debt grid
end

if nargin<4|isempty(dupper)
dupper = 2*(1+rstar); %upper bound of debt grid
end

if nargin<5|isempty(xiMX)
xiMX = 0.99; %elasticity of subsitution between importable and exportables
end

%debt grid
dgrid = linspace(dlower,dupper,nd)';
dd = dgrid(2) - dgrid(1); %debt grid step

n = ny*nd; %total number of states

%Let pX, yN, and d  be 3 ny-by-nd matrices containing values  of the current state
pX = repmat(pXgridlevel,1,nd);  
yN = repmat(yNgridlevel,1,nd);  
d = repmat(dgrid',ny,1);
dix = repmat(1:nd,ny,1); %ny-by-nd matrix of indices of the current debt state
yX = yXendowmentconstant*ones(size(yN));
yM = thet*ones(size(yN));
pX800 = repmat(pX(:),1,800);

%now we just need to adjust these according to the pX so that they have the
%optimal ratio to eachother. Let Z be the optimalratio
changeindebt = bsxfun(@plus,dgrid'/(1+rstar),-d(:));

Z = ((1-b)/b)^xiMX * pX(:).^-xiMX; 

%ADJnoswappingcXtry = bsxfun(@times, noswappingcXtry, repmat(pX(:), 1, 800));

%aux1 = ADJnoswappingcXtry + noswappingcMtry;
%cMtry = bsxfun(@times, aux1, (Z+1).^-1);
MXaux = changeindebt + repmat(yM(:), 1, 800) + bsxfun(@times, repmat(yX(:), 1, 800), pX(:));
testcMtry = MXaux ./ (1 + (pX800).^(1-xiMX))*(((1-b)/b)^xiMX);


%we need to handle what to do when a ex/im constraint would be binding, in this
%case they would consume all exportables or import no additional
%importables 
testcMtry2 = testcMtry;

testcMtry2(testcMtry <= thet) = testcMtry(testcMtry <= thet) - thet;
testcMtry2(testcMtry > thet) = 0;


testcXtry = MXaux ./ ((pX800) + (pX800).^xiMX*((((1-b)/b)^-xiMX)));
testcXtry2 = testcXtry;
testcXtry2(testcXtry >= yXendowmentconstant) = testcXtry(testcXtry >= yXendowmentconstant) - yXendowmentconstant;
testcXtry2(testcXtry < yXendowmentconstant) = 0;
%the optimal ratio only holds when the constraints on importable and
%exportable goods aren't binding, so we need to subtract the difference,
%then convert to the value.
testcMtry = testcMtry + testcXtry2.*(pX800);
testcXtry = testcXtry + (testcMtry2.*(pX800).^-1);

testcMtry(testcMtry < thet) = thet;
testcXtry(testcXtry > yXendowmentconstant) = yXendowmentconstant;

cMtry = testcMtry;
cXtry= testcXtry;

cMtry(cMtry<=0) = NaN;
cXtry(cXtry<=0) = NaN;



%cXtry = bsxfun(@times, cMtry, Z);

clear noswappingcMtry noswappingcXtry necessaryexportables necessaryimportables necessaryexportablesaux

cTtry = ((b*((cMtry).^(1-1/xiMX)))+(1-b)*((cXtry).^(1-1/xiMX))).^(1/(1-1/xiMX));


%now that we have all the cT values, we can use a lot of the code as before
%n-by-nd matrix of relative price of nontradables as a function of the current state (n in total) and possible values for next-period debt (nd in total)
bsxfun(@rdivide, cTtry, yN(:)); 
pNtry  = (1-a)/a*ans.^(1/xi);



%now we have the the prices of nontradables as a function of the states
%we need the collateral constraint now




%the composite consumption good of (M, X, N combined).
bsxfun(@plus, a * cTtry.^(1-1/xi),(1-a) * yN(:).^(1-1/xi));
ctry = (ans).^(1/(1-1/xi));

%marginal utility of tradable consumption, doesn't really make much sense
%but can fix later if needed...
latry = a*ctry.^(-sigg) .* (cTtry./ctry).^(-1/xi);

%initializations
MXasdf = (1-(1/xiMX));
%test = (b*(yXendowmentconstant*min(pX(:))).^MXasdf + (1-b)*(yM(:).^MXasdf)).^(1/(MXasdf));
cT0 = min((b*(yXendowmentconstant*pX(:)).^MXasdf + (1-b)*(yM(:).^MXasdf)).^(1/(MXasdf)))- (rstar / (1+rstar)*d);
c0 = (a * cT0.^(1-1/xi) + (1-a) * yN.^(1-1/xi)).^(1/(1-1/xi));
la0 = a*c0.^(-sigg) .* (cT0./c0).^(-1/xi);
la = la0; %marginal utilliy of tradable consumption
laA = zeros(ny,nd); %marginal utility of tradable consumption when unconstrained

%bsxfun(@times, pNtry, yN(:));
collateraltry = kapa*bsxfun(@plus,((yXendowmentconstant*pX(:))+yM(:)),bsxfun(@times, pNtry, yN(:)));


%subtract the debt to see where it is slack
slacktry = bsxfun(@minus,collateraltry,dgrid');


%find the points we are looking for where there is crisis
aux1 = sign(slacktry);
aux2 = sign([slacktry(:,1) slacktry(:,1:end-1)]);
aux3 = sign([ slacktry(:,2:end) slacktry(:,end)]);
bindCix = find(aux1.*aux2==-1 & slacktry<0); %Points of type C in the graphs of the paper cited in the preamble of this program
bindBix = find(aux1.*aux3==-1 & slacktry<0); %Points of type B in the graphs of the paper cited in the preamble of this program
clear aux1 aux2 aux3
slacktry(bindCix) = 0;
slacktry(bindBix) = 0;

dist = inf;

while  dist>0

%This loop computes  pai*la but allowing for the product pai(i,j)*la(j,k)=0 if pai(i,j)=0 and la(j,k)=NaN
aux = zeros(ny,nd);
for i=1:ny
laaux = la;
pai1 = pai(i,:);
find(pai1==0);
laaux(ans,:) = 0;
aux(i,1:nd) = pai1*laaux;
end

%n-by-nd matrix containing values of  mu_t*la_t as a function of the current state (n in total) and possible values of next-period debt (nd in total)


MUtry = latry/(1+rstar) -  betta*repmat(aux,nd,1);

%Gridize iMUtry=0, by identifying sign change along rows and setting them to 0
aux1 = sign(MUtry);
aux2 = sign([MUtry(:,1) MUtry(:,1:end-1)]);
%MUtry(aux1~=aux2 & ~isnan(aux1.*aux2)) = 0;
%alternative, make the last positive MU already 0
MU1ix = find(aux1.*aux2==-1 & MUtry<0);
clear aux1 aux2 
MUtry(MU1ix) = 0;
clear MU1ix

%Find policy function assuming a nonbinding  collateral constraint
[auxA,dpixA] = min(abs(MUtry),[],2);
pickerA = sub2ind([n nd],(1:n)',dpixA(:));  %In general pickerx (in this case x=u) picks the policy choice out of any n-by-nd matrix  containing for each state (n in total) possible policy choices (nd in total)
laA(:) = latry(pickerA); %marginal utility of tradabels when the collateral constraint is  nonbinding 
find(slacktry(pickerA)<0 | isnan(MUtry(pickerA))  );
laA(ans) = NaN; %rule out violations of the collateral constraint and  negative consumption 
pickerAok = find(~isnan(laA));
muA=MUtry(pickerA);
muA(isnan(laA))=NaN; 

%marginal utility of tradables under the assumption that 'c' points are chosen. 
laC = NaN(ny,nd);
find(MUtry>=0 & ~isnan(MUtry));
pickerCok = intersect(ans,bindCix);
[pickerCiok,~] = ind2sub([n nd],pickerCok);
laC(pickerCiok) = latry(pickerCok);

%marginal utility of tradables under the assumption that 'b' points are chosen. 
laB = NaN(ny,nd);
find(MUtry>=0 & ~isnan(MUtry));
pickerBok = intersect(ans,bindBix);
[pickerBiok,ans] = ind2sub([n nd],pickerBok);
laB(pickerBiok) = latry(pickerBok);

laold = la;

%Construct la under the different equilibrium selection criteria

if strcmp(eqm_selection_criterion,'ab')
la = laC;
la(pickerBiok) = laB(pickerBiok);
la(pickerAok) = laA(pickerAok);
end

if strcmp(eqm_selection_criterion,'ac')
la = laB;
la(pickerCiok) = laC(pickerCiok);
la(pickerAok) = laA(pickerAok);
end

if strcmp(eqm_selection_criterion,'b')
la = laA;
la(pickerCiok) = laC(pickerCiok);
la(pickerBiok) = laB(pickerBiok);
end

if strcmp(eqm_selection_criterion,'c')
la = laA;
la(pickerBiok) = laB(pickerBiok);
la(pickerCiok) = laC(pickerCiok);
end

%Compute the distance between the old and new values of la to determine convergence. This is longer than just computing max(abs(la(:)-laold(:))) because NaN enries must be accounted for

find(isnan(la));
aux1 = la; 
aux1(ans)=0;
find(isnan(laold));
aux2 = laold; 
aux2(ans)=0;
dist = max(abs(aux1(:)-aux2(:)));
end  %while  dist>0









%debt policy function in index form (entries take values between 1 and nd)
dpix  = zeros(ny,nd); 
abs(bsxfun(@minus,latry,la(:)));
clear latry
[~,dpix(:)] = min(ans,[],2);
dpix(isnan(la)) = NaN;

picker = sub2ind([n nd],(1:n)',dpix(:));
pickerok = find(~isnan(picker));

%policy function for the slackness of the collateral constraint
slack = nan(ny,nd);
slack(pickerok) = slacktry(picker(pickerok));
clear slacktry

%equilibrium mu
mu = nan(ny,nd);
mu(pickerok) = MUtry(picker(pickerok));
clear MUtry

%Debt policy function
dp = nan(ny,nd);
dp(pickerok) = dgrid(dpix(pickerok));

%Policy functions for cT, c,  and p
%cT = pX + dp/(1+rstar) -d;

cM = (pX*yXendowmentconstant + yM + dp/(1+rstar) -d) .* (((1-b)/b)^xiMX*(pX.^-xiMX).*pX + 1).^-1;
cX = ((1-b)/b)^xiMX*(pX.^-xiMX).*cM;
cT =  b * cM.^(1-1/xiMX)+ (1-b)*cX.^(1-1/xiMX);
cT = cT.^(1/(1-1/xiMX));
c =  a * cT.^(1-1/xi)+ (1-a)*yN.^(1-1/xi);
c = c.^(1/(1-1/xi));
p = (1-a)/a *(cT./yN).^(1/xi);

%crisis in the present period
crisis = slack;
crisis(slack==0) = 1;
crisis(slack~=0) = 0;
crisis(isnan(slack)) = nan;

%slope of the RHS of the collateral constraint NEED FIX
slope = kapa /xi ./(1+rstar) .* (p./cT) .* yN;


eval(['save ' filename '.mat filename nd ny n betta sigg bindBix bindCix   pai xi a dgrid pXgrid yNgrid pXgridlevel yNgridlevel  laB laC laA  la dpix dp  crisis pX yN cM cX cT c kapa  picker  dupper dlower dd rstar   d p  slack mu  bindCix bindBix   slope  '])
