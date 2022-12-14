function data = sim_pan_probit(n,D,mod,x);
% simulate panel probit data using model mod 
% for n individuals and D decisions for each.
%
% INPUT: n ... integer, number of inidividuals.
%        D ... integer, number of choice dec per indiv.
%        mod ... pan_prob_mod.
%        [x] ... regressor array.
%
% OUTPUT data ... vector of choice data.
%
% AUTHOR: dbauer, 24.9.2015. 

% get parameters from the model. 
b = mod.b;
O = mod.O;
L = mod.L;

K = mod.k;
M = mod.m; 

if nargin<4
    x = randn(n,D,M,K);
end;

opt = ones(D,M);

for i=1:n % cycle over individuals.
    bi = b + O* randn(K,1); % individ. specific coefficients. 
    for d=1:D % cycle over decisions. 
        Uid = squeeze(x(i,d,:,:))*bi + L*randn(M,1); % vector of random utils. 
        [~,y(d)] = max(Uid); % decision maximizes utility. 
    end;
    
    data(i)= choice_data(y,squeeze(x(i,:,:,:)),opt);
end; 
        