function newval = cdfmvna_SJ2(a,r,sys_randper);
% cdfmvna Solow Joe method.
% dbauer, 29.6.2015.

if nargin<3
    sys_randper = 1;
end;

persistent w 

m = length(a); % dimension of CDF.

%s1 = 0; % w...permutation of 1:m, s1: new seed, in case approx is negative.
if(sys_randper ==0)
    if isempty(w)
        w = getComb(m,2); % all combinations of two dims ordered, rest unordered.
    end;
%    s1 = s;
else
    w = zeros(sys_randper,m);
    rng(1232);
    for z =1:sys_randper
        w(z,:) = randperm(m);
    end
end;
n1 = size(w,1); % num permutations.
lp = 0; % perform all calcs in logs.

%	a <- 	a * (a < 5.7) + 5.7*(a >= 5.7); % truncate a at 5.7.
%	ab = t(a);
for j =1:n1
    x = a(w(j,:));
    rho = r(w(j,:),w(j,:));
    rho(rho>1)=1;
    rho(rho<-1)=-1;
    z=zeros(m,m);
    for k=1:m
        for l=(k+1):m
            %z(k,l)= mvnxpb(rho([k,l],[k,l]),-Inf*[1;1],[x(k),x(l)]');
            z(k,l)= bvnu( -x(k), -x(l), rho(k,l)); 
            z(l,k)=z(k,l);
        end;
    end;
    
    z3 = normcdf(x); % univariate.
    z = z+ diag(z3); % replace bivariate by univariate on main diagonal.
    
   
    %ncdfx = normcdf(x);
    z1 = z3* z3'; %ncdfx*ncdfx'; % cdf under independence.
    z2 = 1-z3; % alt. prob. univariate.
    
    lcond =  log(max(z(1,2),0.000000000000001));
    
    % direct calc of inverse in order to avoid probs with nonsingularity.
   Ainv = zeros(m,m);
    M = z(1:2,1:2) - z1(1:2,1:2);
    detM = max(det(M),10^(-6));
    Ainv(1:2,1:2) = [M(2,2),-M(1,2);-M(2,1),M(1,1)]/detM;

    for k = 3:m % here is the approximation!!
        omega21 = z(k,1:k-1)-z1(k,1:k-1); % equ. (7)
        lcondk = log(max(0.000000000000001,z3(k) + omega21 *(Ainv(1:(k-1),1:(k-1)))*z2(1:k-1))) ;   % equ (8).
        lcond=lcond+lcondk;
        % update
        B = Ainv(1:(k-1),1:(k-1)) * omega21';
        update = [-B;1];
        Ainv(1:k,1:k) = Ainv(1:k,1:k) + (update * update' )/(z(k,k)-z1(k,k)-omega21*B);
    end;
    
    %lpcomb = lz+lcond; % new estimate.
    lp = lp+lcond;  % add to calculate average.
end;
newval = [lp/n1];

        