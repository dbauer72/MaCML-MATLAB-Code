function [newval,grad] = mcdfmvna_SJ2(a,r,sys_randper)
% cdfmvna Solow Joe method.
% log(prob) and gradient.
% dbauer, 23.9.2015.

m = length(a); % dimension of CDF.
w = 0;

if nargin<3
%     if m<5 
%         sys_randper = 0;
%     else
%         sys_randper = 1;
%     end;
    sys_randper = 1;
end;


%s1 = 0; % w...permutation of 1:m, s1: new seed, in case approx is negative.
if(sys_randper ==0)
    w = getComb(m,2); % all combinations of two dims ordered, rest unordered.
    %    s1 = s;
else
    % draw randperm once for ever.
    w = 1:m;
%     w = zeros(sys_randper,m);
%     rng('shuffle');
%     for z =1:sys_randper
%         w(z,:) = randperm(m);
%     end
end
n1 = size(w,1); % num permutations.
lp = 0; % perform all calcs in logs.

gw = zeros(1,m);
gr = zeros(1,m*(m-1)/2);

[I,J]= find(tril(ones(m),-1));

for j=1:n1 % cycle over permutations.
    % initialize vars.
    x = a(w(j,:));
    rho = r(w(j,:),w(j,:));
    rho(rho>1)=1;
    rho(rho<-1)=-1;
    z=zeros(m,m);
    for k=1:m
        for l=(k+1):m
            z(k,l)= bvnu( -x(k), -x(l), rho(k,l) );
            %mvnxpb_new(rho([k,l],[k,l]),[x(k),x(l)]'); % evaluate only at upper bound, diagonal equal to one. 
            z(l,k)=z(k,l);
        end;
    end;
    
    Phi = phid(x); % univariate.
    phi = phip(x);
    z = z+ diag(Phi); % replace bivariate by univariate on main diagonal.
    
    %ncdfx = phid(x);
    z1 = Phi* Phi'; %ncdfx*ncdfx'; % cdf under independence.
    X0 = 1-Phi; % alt. prob. univariate.
    lcond =  log(max(z(1,2),0.000000001));
    
    % direct calc of inverse in order to avoid probs with nonsingularity.
    Ainv = zeros(m,m);
    M = z(1:2,1:2) - z1(1:2,1:2);
    detM = max(det(M),10^(-6));
    Ainv(1:2,1:2) = [M(2,2),-M(1,2);-M(2,1),M(1,1)]/detM;
    
    % initialize for gradient.
    gw1 = zeros(1,m);
    gr1 = zeros(1,m*(m-1)/2);
    
    phi = phip(x); % PDF at x.
    muygx = rho*0;
    
    rho1 = rho;
    rho1 = rho1- diag(diag(rho1));
    rho2 = sqrt(1-rho1.^2);
    rho2 = rho2- diag(diag(rho2))+eye(m);
    muygx = (ones(m,1)*x' - rho.*(x*ones(1,m)))./rho2;
    nmuygx = phid(muygx);
    dPhixy = diag(phi)*nmuygx; % pdf x * prob of stud. resid.
    phiPhi = phi * Phi'; % pdf x * prob of y.
    dCov = dPhixy-phiPhi; % derivative of Cov(I_i,I_j) w.r.t. w_i.
    dCov = dCov - diag(diag(dCov)) + diag(phi.*(1-2*Phi));
    
    dPhi2 = (diag(phi)*phip(muygx))./rho2;
    dPhi2 = dPhi2- diag(diag(dPhi2));
    
    for l=1:length(I)
        dPhi2Vec(l) =  dPhi2(I(l),J(l));
    end;
    
    kk=length(I); % number of corrs.
    
    % gradient calcs.
    gw1(2) = dPhixy(2,1)/exp(lcond);
    gw1(1) = dPhixy(1,2)/exp(lcond);
    gr1(1) = dPhi2Vec(1)/exp(lcond);
    
    for k=3:m  % calc of approximate CDF vals.
        omega21 = z(k,1:k-1)-z1(k,1:k-1); % equ. (7
        invomg11 = Ainv(1:(k-1),1:(k-1));
        ioX0 = (invomg11 * X0(1:k-1));
        oio11 = (omega21 * invomg11);
        condk = Phi(k) + omega21 *ioX0;
        lcondk = log(max(0.000000000001,condk));   % equ (8).
        lcond=lcond+lcondk;
        
        
        %%%   Start here for gradients with respect to rho parameters   ###
        gw2 = zeros(1,m);
        gr2 = zeros(1,kk);
        
        % deriv w.r.t. w_i's., i<k
        
%         gw2(1:(k-1)) = (ioX0 .* dCov(1:(k-1),k))';
%         for ii=1:(k-1)
%             gw2(ii)=gw2(ii) - (dCov(ii,1:(k-1))* ioX0) *oio11(ii) - (oio11* dCov(ii,1:(k-1))')*ioX0(ii) + ioX0(ii)*oio11(ii)*dCov(ii,ii);
%         end
%         gw2(1:(k-1))= gw2(1:(k-1)) - (oio11(1:(k-1)) .* phi(1:(k-1))');
        gw2(1:(k-1)) = (ioX0 .* dCov(1:(k-1),k)) - oio11(:).*(dCov(1:(k-1),1:(k-1))*ioX0) - ioX0.*(dCov(1:(k-1),1:(k-1))*oio11(:)) + ioX0 .* oio11(:) .* diag(dCov(1:(k-1),1:(k-1)))- (oio11(1:(k-1))' .* phi(1:(k-1)));
        
%         if norm(gw2(1:k-1)-gw3')>10^(-6)
%             disp('here')
%         end;
        
        % deriv w.r.t w_k.
        gw2(k)=phi(k)+ dCov(k,1:(k-1)) * ioX0;
        
        % deriv w.r.t. rho's.
        cur = 1;
        for ii = 1:(k-1)
            if (ii<k-1)
                for jj = (ii+1):(k-1)
                    gr2(cur) = gr2(cur) + (- ioX0(ii)*oio11(jj)-ioX0(jj)*oio11(ii))*dPhi2(jj,ii);
                    cur=cur+1;
                end
            end
            gr2(cur)=gr2(cur)+dPhi2(k,ii)*ioX0(ii);
            cur = cur + 1 + (m-k);
        end
        
        % add to overall function.
        gr2  = gr2/condk;
        gw2 = gw2/condk;
        
        gr1 = gr1+gr2;
        gw1 = gw1+gw2;
        
        % update
        B = Ainv(1:(k-1),1:(k-1)) * omega21';
        update = [-B;1];
        Ainv(1:k,1:k) = Ainv(1:k,1:k) + (update * update' )/(z(k,k)-z1(k,k)-omega21*B);
    end % k=3:m-1.
    
    %    commands below to resequence gradients based on permutation   #####
    
    [~,ri]= sort(w(j,:));

    gw1(w(j,:))=gw1;

    rr = 0*rho;
    for l=1:length(I)
        rr(I(l),J(l))=l;
        rr(J(l),I(l))=l;
    end;
    rr = rr(ri,ri);
    
    for l=1:length(I)
        Ir(l)=rr(I(l),J(l));
    end;
    
    gr1=gr1(Ir);
    gw = gw+gw1;
    gr= gr+gr1;
    
    lp = lp+lcond;  % sum of all prob values. is only provided with output, no other use.
end % cycle over permutations.

newval = [lp/n1];
grad = [gw,gr]/n1;




    


