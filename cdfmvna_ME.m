function [lp,dlp] = cdfmvna_ME(Zj,R); 
% cdfmvna_ME: evaluate approximate log(CDF) according to Mendell Elston.
% a: points where CDF is evaluated.
% r KxK correlation matrix. 
%
% dbauer, 22.9.2015.

cutoff=6;

Zj(Zj>cutoff)=cutoff;
Zj(Zj<-cutoff)=-cutoff;

persistent Phi phi

if isempty(Phi)
    xgr = [-6:0.0001:6];
    Phi = phid(xgr);
    phi = phip(xgr);
    clear xgr;
end;

K = length(Zj); % dimension of CDF.
% reorder 
[Zj,i]=sort(Zj,'descend');

Rij = R(i,i);
P = Phi(min(max(1,floor((Zj(1)+6)*10000)),length(Phi)));  % 
p = phi(min(max(1,floor((Zj(1)+6)*10000)),length(phi))); %
lp = log(P); % perform all calcs in logs.

if nargout>1 % gradient also to be calculated 
    % initialize. 
    [I,J]= find(tril(ones(K),-1));
    npar = K+length(I); % first K entries deriv w.r.t. Zj, rest w.r.t. entries in R. 
    dRij = zeros(K,K,npar);
    dtRij = dRij;
    for j=1:length(I)
        dRij(I(j),J(j),K+j)=1;
        dRij(J(j),I(j),K+j)=1;
    end;
    
    dZj = [eye(K),zeros(K,npar-K)];
    dlp = dZj(1,:)*p/P;
    dtZ = dZj*0;
    
    dsig = dZj*0;
    dajjm1 = dlp*0; 
end;
for j = 1:(K-1)
    ajjm1 = p/P;
    tZ = Zj + ajjm1*Rij(:,j);    % update Zj  
    R_jj = (Rij(:,j)*Rij(j,:));
    tRij = Rij - R_jj*(ajjm1+Zj(j)) *ajjm1;    % update Rij
    sig = sqrt(diag(tRij)); % variances
    sig = max(sig,0.00000001);
    if nargout>1
        dajjm1 = -dZj(j,:)*(ajjm1*Zj(j)+ajjm1^2);
        dtZ = dZj + Rij(:,j)*dajjm1 + reshape(dRij(:,j,:),[K,npar])*ajjm1;
        ajjm1p =(ajjm1+Zj(j))*ajjm1;
        R_jj_fac = R_jj * (2*ajjm1 + Zj(j));

        dtRij =  dRij - ajjm1p * (bsxfun(@times,dRij(:,j,:),Rij(j,:))+ bsxfun(@times,Rij(:,j),dRij(j,:,:))) - bsxfun(@times,R_jj_fac,reshape(dajjm1,[1 1 npar])) - ajjm1 * bsxfun(@times,R_jj,reshape(dZj(j,:),[1 1 npar]));
        % 
        for l=1:K,
            dsig(l,:) = -0.5*(dtRij(l,l,:))/(sig(l)^3);
        end;
    end;

    
    Zj = tZ./sig;
    Rij = tRij./(sig*sig'); 
    Rij = Rij-diag(diag(Rij))+eye(K);
    ds = 1./sig;
    DS = diag(ds);
    
    if nargout>1
        DSS = (ds*ds');
        dZj = DS* dtZ + diag(tZ)*dsig;
        % alternative calc. 
        dsdsig(1,:,:)= dsig;
        ppp = bsxfun(@times,ds,dsdsig);
        dRij2 = bsxfun(@times,dtRij,DSS) + bsxfun(@times,tRij,ppp) + bsxfun(@times,tRij,permute(ppp,[2 1 3]));

    end;
    
    P = Phi(min(max(1,floor((Zj(j+1)+6)*10000)),length(Phi))); 
    p = phi(min(max(1,floor((Zj(j+1)+6)*10000)),length(phi))); 
   
    lp = lp + log(P);
    
    if (nargout>1) % gradient calc?
        dlp = dlp + (p/P * dZj(j+1,:));
    end;
end;


if (nargout>1)
    % reorder according to entries.
    [~,ri]= sort(i(:));
    dlp(1:K)=dlp(ri);
    rr = 0*R;
    for l=1:length(I)
        rr(I(l),J(l))=l;
        rr(J(l),I(l))=l;
    end;
    rr = rr(ri,ri);
    for l=1:length(I)
        Ir(l)=rr(I(l),J(l));
    end;
    
    dlp((K+1):end)=dlp(K+Ir);
end;

