function [lp,gr,lpi] = cal_prob_like_cml(mod,data,fun,theta)
% calculates the likelihood and its gradient.
% INPUT:  mod ... pan_prob_mod
%         data ... vector of choice data
%         fun ... string specifying P(.) function incl. derivatives.
%         [theta ... parameter vector]
%
% OUTPUT: lp ... log of likelihood.
%         gr ... gradient of log-like.
%
% AUTHOR: dbauer, 24.9.2015.

if nargin ==4
    mod = pan_prob_mod(theta,mod);
    if nargout>1
        [lp,gr,lpi] = cal_prob_like_cml(mod,data,fun);
    else
        lp = cal_prob_like_cml(mod,data,fun);
    end;
    return
end;
% 


% integers specifying the setting.
n = length(data); % number of inidividuals.
M = mod.m;
K = mod.k;
nvarb = mod.nvarb;
nvarO = mod.nvarO;

Om = mod.Omega;
dOmega = mod.dOmega;
MSigma = mod.MSigma;
dMSigma = mod.dMSigma;
db = mod.db;

b = mod.b;

pars = mod.nvarO+mod.nvarL;
lp = 0;
gr = zeros(1,nvarb + pars);

for i=1:n
    % get data
    D = data(i).D;
    yd = data(i).y;
    dxd = data(i).dx;
    opt = data(i).opt;
    
    nO = sum(opt,2)-1;
    if std(nO)==0
        co = 0;
    else
        co = 1;
    end;
    
    gri = gr*0;
    lpi(i)=0;
    
    if D>1
        comb = nchoosek(1:D,2); % all possible combinations of D decisions.
        for c= 1:size(comb,1)
            % calc kk1, kk2.
            T = nO(comb(c,1))+nO(comb(c,2)); % number of options avail - numbre of dec.
            kk1 = zeros(T,1);
            kk2 = zeros(T,T);
            
            dX = kk1*ones(1,K);
            d_Sig = zeros(T,T);
            gr_Sig = zeros(T,T,pars);
            
            cur = 0;
            for d=comb(c,:) % cycle over indiv. decisions
                if co
                    in = find(opt(d,setdiff(1:M,yd(d))));
                else
                    in = 1:M-1;
                end;
                dX(cur+[1:nO(d)],:)=reshape(dxd(d,in,:),nO(d),size(dX,2));
                d_Sig(cur+[1:nO(d)],cur+[1:nO(d)])=reshape(MSigma(in,in,yd(d)),nO(d),nO(d)); % squeeze(MSigma(in,in,yd(d)));
                gr_Sig(cur+[1:nO(d)],cur+[1:nO(d)],nvarO+1:end)=dMSigma(in,in,yd(d),:);
                cur = cur+length(in);
            end
            
            kk1 = -dX*b;
            kk2 = d_Sig + dX*Om*dX';
            sd = sqrt(diag(kk2));
            
            for o=1:mod.nvarO
                gr_Sig(:,:,o)=dX*reshape(dOmega(:,:,o),K,K)*dX';
            end;
            
            % normalize kk1 and kk2.
            [dkk1,dkk2,nkk2] = calc_deriv_kk2(kk1,kk2,sd,gr_Sig);
            
            % calc probs
            [lpii,grii] = feval(fun,kk1./sd,nkk2);
            lpi(i) = lpi(i) + lpii;
            
            % calc deriv
            gri(1:nvarb)= gri(1:nvarb) - (grii(1:T)./sd')*dX*db; % der. w.r.t. entries in b.
            gri(nvarb+1:end)= gri(nvarb+1:end) +grii(1:T)*dkk1 + grii(T+1:end)*dkk2; % der. w.r.t. entries in O and L.
        end % cycle over different combs.
        lpi(i)=lpi(i)/size(comb,1);
        gri=gri/size(comb,1);
    else % only one decision.
        
        % calc kk1, kk2.
        nO = sum(opt)-1;
        if nO==M-1
            co = 0;
        else
            co = 1;
        end;
        T = nO; % number of options avail - numbre of dec.
        kk1 = zeros(T,1);
        kk2 = zeros(T,T);
        
        dX = kk1*ones(1,K);
        d_Sig = zeros(T,T);
        gr_Sig = zeros(T,T,pars);
        
        
        if co
            in = find(opt(d,setdiff(1:M,yd)));
        else
            in = setdiff(1:M,yd);
        end;
        dX([1:nO],:)=squeeze(dxd(1,in,:));
        d_Sig([1:nO],cur+[1:nO])=squeeze(MSigma(in,in,yd));
        gr_Sig(cur+[1:nO],cur+[1:nO(d)],nvarO+1:end)=dMSigma(in,in,yd,:);
        
        
        kk1 = -dX*b;
        kk2 = d_Sig + dX*Om*dX';
        sd = sqrt(diag(kk2));
        
        for o=1:mod.nvarO
            gr_Sig(:,:,o)=dX*squeeze(dOmega(:,:,o))*dX';
        end;
        
        % normalize kk1 and kk2.
        [dkk1,dkk2,nkk2] = calc_deriv_kk2(kk1,kk2,sd,gr_Sig);
        
        % calc probs
        [lpii,grii] = feval(fun,kk1./sd,nkk2);
        lpi(i) = lpi(i) + lpii;
        
        % calc deriv
        gri(1:nvarb)= gr(1:nvarb) - (gri(1:T)./sd')*dX*db; % der. w.r.t. entries in b.
        gri(nvarb+1:end)= gr(nvarb+1:end) +grii(1:T)*dkk1 + grii(T+1:end)*dkk2; % der. w.r.t. entries in O and L.
        
    end % only one decision.
    if lpi(i)==-Inf
        lpi(i)=-50;
    end;
    lp = lp + lpi(i);
    if ~isnan(gri)
        gr = gr + gri;
    end;
end;
lp = -lp/n;
gr = -gr/n;


function [dkk1,dkk2,nkk2] = calc_deriv_kk2(kk1,kk2,sd,gr_Sig);
% calculates the derivative of normalized kk2 based on the original. 

nkk2 = kk2./(sd*sd');

% this part only works for constant size of matrices!!!
persistent Ir  Jr out

M = size(kk2,1); 
npar = size(gr_Sig,3); 

if isempty(Ir) 
    
    [Ir,Jr]= find(tril(ones(M),0));
    in = 1:M^2;
    for j=1:length(Ir)
        out = [out,(Ir(j)-1)*M+Jr(j)];
    end;
end;

dkk2 = zeros(M^2,npar);
dkk1 = zeros(M,npar);

dind = [1:(M+1):M^2];
vgr_sig = reshape(gr_Sig,[M^2,npar]);
dkk1 = -bsxfun(@times,vgr_sig(dind,:),(kk1./(sd.^3)))/2;

dkk2l_h2 = bsxfun(@times,1./(sd*sd'),gr_Sig);
vdkk2l_h2 = reshape(dkk2l_h2,M^2,npar);
ddkk2l_h2 = vdkk2l_h2(dind,:); 

dkk2 = vdkk2l_h2 - 0.5*bsxfun(@times,reshape(nkk2,M^2,1),ones(1,npar)).*(repmat(ddkk2l_h2,[M,1])+kron(ddkk2l_h2,ones(M,1)));

dkk2(out,:)=[];


