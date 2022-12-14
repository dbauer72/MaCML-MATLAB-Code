% script for testing the procs. 
% dbauer.

n = 500; % number of individuals
D = 5; % number of decisions per individual
M = 5; % number of options per decision
K = 5; % number of regressors. 
IT = 100;


fun_ME = 'cdfmvna_ME';
fun_SJ = 'mcdfmvna_SJ';

% generate model according to Bhat standard model. 
Hb = [zeros(1,4);eye(4)];
fb = [1.5,-1,2,1,-2];

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
fO([1,7,13,19,25])=1;

theta = [theta,thO]; 


HL = zeros(25,0);
fL = zeros(25,1);

mod = pan_prob_mod(theta(1:4)',Hb,fb',theta(5:end)',HO,fO,zeros(1,0),HL,fL,{'O1','O2','O3','O4','O5'},{'R1','R2','R3','R4','R5'});

options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on','Display','iter');
clear thce_ME thce_SJ thce_MXL

% mixed logit estimation
lm2 = cMIXLOGIT();
lm2.B = fb;
lm2.IDV = [1:5;ones(1,5)]'
lm2.COV = ones(5,1);
lm2.IDV = [4:8;ones(1,5)]';
lm2.NAMES = {'O1','O2','O3','O4','O5'}';
lm2.B = fb';
lm2.ASC = [];
Omeg. fO= mod.fO;
Omeg.HO = mod.HO;
Omeg.thO = mod.thO;
lm2.W = Omeg;
lm2.NDRAWS = 200;

h = waitbar(0,'Please wait...'); 

for t=1:IT % generate 100 instances. 
    waitbar(t/IT,h);
    % generate data
    %rng(rem(fix(now*1000000),1000));
    data = sim_pan_probit(n,D,mod);

    % CML estimation with ME
%     tic;
%     f_ME = @(x) cal_prob_like_cml(mod,data,fun_ME,x);
%     thce_ME(t,:) = fminunc(f_ME,theta,options);
%     time_ME(t) = toc; 
%     
    % CML estimation with SJ
    tic;
    f_SJ = @(x) cal_prob_like_ccml(mod,data,fun_SJ,x);
    thcce_SJ(t,:) = fminunc(f_SJ,theta,options);
    time_ccSJ(t) = toc; 
    
    % CCML estimation with ME
    tic;
    f_ME = @(x) cal_prob_like_ccml(mod,data,fun_ME,x);
    thcce_ME(t,:) = fminunc(f_ME,theta,options);
    time_ccME(t) = toc; 
    
    % Mixed Logit estimation. 
%     tic;
%     MLe = fit_mixlogit(data,lm2);
%     thce_MXL(t,:)= MLe.PARAMHAT(:)'; 
%     time_MXL(t) = toc; 
end;

close(h); 