% script for testing the MaCML estimation procedures.
% 
% dbauer, 13.12.2022.

n = 500; % number of individuals
D = 5; % number of decisions per individual
M = 5; % number of options per decision
K = 5; % number of regressors. 

%fun = 'cdfmvna_ME'; % select Solow Joe approximation
 fun = 'mcdfmvna_SJ'; % alternatively select Mendel Elston. 

% generate model according to Bhat standard model. 
Hb = [zeros(1,4);eye(4)];
fb = [1,0,0,0,0];

theta=[-1,2,1,-2];

HO = zeros(25,15);
c =0;
thO = [];
for i=1:5
    for j=i:5
        c = c+1;
        HO((i-1)*5+j,c)=1;
        if (j==i)
            thO(c)=1;
        else
            thO(c)=0;
        end;
    end;
end;

fO=zeros(25,1);
%fO([1,7,13,19,25])=1;

theta = [theta,thO]; 

HL = zeros(25,0);
fL = zeros(25,1);
fL([1,7,13,19,25])=1;

for j=1:5
    Opt_Names{j} = num2str(j);
    Reg_Names{j} = sprintf('X%d',j);
end;

tic;
mod = pan_prob_mod(theta(1:4),Hb,fb',theta(5:end)',HO,fO,zeros(1,0),HL,fL,Opt_Names,Reg_Names);

% generate data
data = sim_pan_probit(n,D,mod);

% estimate model using full pairwise CML
f = @(x) cal_prob_like_cml(mod,data,fun,x);
options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on','Display','iter');

the = fminunc(f,theta,options);
toc
norm(theta-the)

% estimate model using adjacent pairwise CML
tic;
mod = pan_prob_mod(theta(1:4),Hb,fb',theta(5:end)',HO,fO,zeros(1,0),HL,fL,Opt_Names,Reg_Names);

% generate data
data = sim_pan_probit(n,D,mod);

% estimate model using full pairwise CML
f = @(x) cal_prob_like_ccml(mod,data,fun,x);

thce = fminunc(f,theta,options);
toc
norm(theta-thce)



