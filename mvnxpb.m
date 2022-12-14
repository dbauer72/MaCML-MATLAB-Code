function  pb = mvnxpb( r, a, b )
%  p = mvnxpb( r, a, b )
%    uses bivariate conditionally expected values for variables to estimate 
%    an MVN probability for positive definite covariance matrix r,
%    with lower integration limits a and upper integration limits b. 
%    Covariance matrix R is permuted to improve accuracy
%   Probability p is output
%    Example:
%     r = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5];
%     a = -inf(4,1); b = [ 1 2 3 4 ]';
%     p = mvnxpb( r, a, b ); disp(p)
%
%   Alan Genz is the author of this function and following Matlab functions.
%     Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
%     Email : AlanGenz@wsu.edu
%
%
[n,n] = size(r); c = r; ap = a; bp = b; ep = 1e-10; % singularity tolerance;
d = sqrt(max(diag(c),0)); 
% for i = 1 :  n  % scale c to make it a correlation matrix
%   if d(i) > 0, ap(i) = ap(i)/d(i); bp(i) = bp(i)/d(i);
%     c(:,i) = c(:,i)/d(i); c(i,:) = c(i,:)/d(i);     
%   end
% end
% Dynamically determine Cholesky factor, 
%    permuting variables to minimize outer integrals
y = zeros(n,1); p = 1; pb = 1; D = eye(n);
for k = 1 : n, im = k; ckk = 0; dem = 1; s = 0; 
  for i = k : n 
    if c(i,i) > eps, cii = sqrt( max( [c(i,i) 0] ) ); 
      if i > 1, s = c(i,1:k-1)*y(1:k-1); end
      ai = ( ap(i)-s )/cii; bi = ( bp(i)-s )/cii; de = Phi(bi) - Phi(ai);
      if de <= dem, ckk = cii; dem = de; am = ai; bm = bi; im = i; end
    end
  end
  if im > k, km = [1:k-1]; ip = [im+1:n]; ki = [k+1:im-1];
    ap([im k]) = ap([k im]); bp([im k]) = bp([k im]); 
    c([im k],km) = c([k im],km); c(ip,[im k]) = c(ip,[k im]);
    t = c(ki,k); c(ki,k) = c(im,ki)'; c(im,ki) = t'; c(im,im) = c(k,k); 
  end
  if ckk > ep*k, c(k,k) = ckk; c(k,k+1:n) = 0; 
    for i = k+1 : n, c(i,k) = c(i,k)/ckk; 
      c(i,k+1:i) = c(i,k+1:i) - c(i,k)*c(k+1:i,k)';
    end
    if abs(dem) > ep, y(k) = ( phi(am) - phi(bm) )/dem; 
    else, y(k) = ( am + bm )/2;
      if am < -9, y(k) = bm; elseif bm > 9, y(k) = am; end
    end
  else, c(k:n,k) = 0; y(k) = ( ap(k) + bp(k) )/2;
  end, p = p*dem;
  % Compute bivariate product
  if mod(k,2) == 0, u = c(k-1,k-1); v = c(k,k); w = c(k,k-1); 
    c(k-1:n,k-1:k) = c(k-1:n,k-1:k)*[1/u 0;-w/(u*v) 1/v]; 
    ab = ap(k-1:k); bb = bp(k-1:k); cb = 0;
    if k > 2, cb = c(k-1:k,1:k-2)*y(1:k-2); end
    sg = [u^2 u*w;u*w w^2+v^2]; D(k-1:k,k-1:k) = sg;  
    [Ex Ey bv] = bvnmmg( ab-cb, bb-cb, sg ); 
    pb = pb*bv; y(k-1:k) = [ Ex; Ey ]; 
  end
end, if mod(n,2) == 1, pb = pb*dem; end
%
% end mvnxpb
%
function p = Phi(z), p = erfc( -z/sqrt(2) )/2;   % Normal cdf
function p = phi(z), p = exp(-z^2/2)/sqrt(2*pi); % Normal pdf
function [Ex Ey p] = bvnmmg( a, b, sg )
%
%  bvnmmg( a, b, sg )
%  A function for computing bivariate normal probability moments;
%  bvnmmg calculates expected values Ex, Ey and probability p, for 
%    bivariate normal x, with a < x < b and covariance matrix sg.
%  Example: 
%   [Ex Ey p] = bvnmmg([-inf -2],[5 inf],[4 1;1 3]); disp([Ex Ey p])
%
  cx = sqrt(sg(1,1)); cy = sqrt(sg(2,2)); r = sg(2,1)/(cx*cy); 
  xl = a(1)/cx; xu = b(1)/cx; yl = a(2)/cy; yu = b(2)/cy;
  [Ex Ey p] = bvnmom( xl, xu, yl, yu,  r ); Ex = Ex*cx; Ey = Ey*cy;
% end bvnmmg
function [Ex Ey p] = bvnmom( xl, xu, yl, yu,  r ), rs = 1/sqrt(1-r^2);
%
%  A function for computing bivariate normal probability moments;
%  bvnmom calculates expected values Ex, Ey, for bivariate normal x, y,
%   with xl < x < xu and yl < y < yu, and correlation coefficient r.
%  
%  This function uses generalizations of formulas found in   
%    Moments of a Truncated Bivariate Normal Distribution
%    S. Rosenbaum, JRSS B 23 (1961) pp. 405-408.
%
%
p = bvnu(xl,yl,r) - bvnu(xu,yl,r) - bvnu(xl,yu,r) + bvnu(xu,yu,r); 
if xl == -inf && yl == -inf && xu == inf && yu == inf, Ex = 0; Ey = 0;
elseif xl == -inf && xu ==  inf && yl == -inf, Ey = -phi(yu); Ex = 0;
elseif xl == -inf && xu ==  inf && yu ==  inf, Ey =  phi(yl); Ex = 0;
elseif xl == -inf && yl == -inf && yu ==  inf, Ex = -phi(xu); Ey = 0;
elseif xu ==  inf && yl == -inf && yu ==  inf, Ex =  phi(xl); Ey = 0;
elseif xl == -inf && xu ==  inf,     Ey = phi(yl) - phi(yu); Ex = 0;
elseif yl == -inf && yu ==  inf,     Ex = phi(xl) - phi(xu); Ey = 0;
else
  if xl == -inf && yl == -inf 
    Phiyxu = Phi( (yu-r*xu)*rs); pPhixy = -phi(xu)*Phiyxu; 
    Phixyu = Phi( (xu-r*yu)*rs); pPhiyx = -phi(yu)*Phixyu; 
  elseif xu ==  inf && yu ==  inf
    Phiyxl = Phi(-(yl-r*xl)*rs); pPhixy =  phi(xl)*Phiyxl; 
    Phixyl = Phi(-(xl-r*yl)*rs); pPhiyx =  phi(yl)*Phixyl; 
  elseif xl == -inf && yu ==  inf
    Phiyxu = Phi(-(yl-r*xu)*rs); pPhixy = -phi(xu)*Phiyxu; 
    Phixyl = Phi( (xu-r*yl)*rs); pPhiyx =  phi(yl)*Phixyl; 
  elseif xu ==  inf && yl == -inf
    Phiyxl = Phi( (yu-r*xl)*rs); pPhixy =  phi(xl)*Phiyxl;
    Phixyu = Phi(-(xl-r*yu)*rs); pPhiyx = -phi(yu)*Phixyu;
  elseif xl == -inf
    Phiyxu = Phi((yu-r*xu)*rs) - Phi((yl-r*xu)*rs); pPhixy = -phi(xu)*Phiyxu;
    Phixyl = Phi( (xu-r*yl)*rs); Phixyu = Phi( (xu-r*yu)*rs);
    pPhiyx = phi(yl)*Phixyl - phi(yu)*Phixyu; 
  elseif xu ==  inf
    Phiyxl = Phi((yu-r*xl)*rs) - Phi((yl-r*xl)*rs); pPhixy =  phi(xl)*Phiyxl;
    Phixyl = Phi(-(xl-r*yl)*rs); Phixyu = Phi(-(xl-r*yu)*rs);
    pPhiyx = phi(yl)*Phixyl - phi(yu)*Phixyu; 
  elseif yl == -inf
    Phiyxl = Phi( (yu-r*xl)*rs); Phiyxu = Phi((yu-r*xu)*rs);
    pPhixy = phi(xl)*Phiyxl - phi(xu)*Phiyxu; 
    Phixyu = Phi((xu-r*yu)*rs) - Phi((xl-r*yu)*rs); pPhiyx = -phi(yu)*Phixyu;
  elseif yu ==  inf
    Phiyxl = Phi(-(yl-r*xl)*rs); Phiyxu = Phi(-(yl-r*xu)*rs);
    pPhixy = phi(xl)*Phiyxl - phi(xu)*Phiyxu; 
    Phixyl = Phi((xu-r*yl)*rs) - Phi((xl-r*yl)*rs); pPhiyx = phi(yl)*Phixyl;
  else  
    Phiyxl = Phi((yu-r*xl)*rs) - Phi((yl-r*xl)*rs);
    Phiyxu = Phi((yu-r*xu)*rs) - Phi((yl-r*xu)*rs);
    pPhixy = phi(xl)*Phiyxl - phi(xu)*Phiyxu;
    Phixyl = Phi((xu-r*yl)*rs) - Phi((xl-r*yl)*rs);
    Phixyu = Phi((xu-r*yu)*rs) - Phi((xl-r*yu)*rs);
    pPhiyx = phi(yl)*Phixyl - phi(yu)*Phixyu;
  end, Ex = pPhixy + r*pPhiyx; Ey = r*pPhixy + pPhiyx;
end, Ex = Ex/p; Ey = Ey/p;
%
%  end bvnmom
%
function p = bvnu( dh, dk, r )
%BVNU
%  A function for computing bivariate normal probabilities.
%  bvnu calculates the probability that x > dh and y > dk. 
%    parameters  
%      dh 1st lower integration limit
%      dk 2nd lower integration limit
%      r   correlation coefficient
%  Example: p = bvnu( -3, -1, .35 )
%  Note: to compute the probability that x < dh and y < dk, 
%        use bvnu( -dh, -dk, r ). 
%

%   Author
%       Alan Genz
%       Department of Mathematics
%       Washington State University
%       Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%
%    This function is based on the method described by 
%        Drezner, Z and G.O. Wesolowsky, (1989),
%        On the computation of the bivariate normal inegral,
%        Journal of Statist. Comput. Simul. 35, pp. 101-107,
%    with major modifications for double precision, for |r| close to 1,
%    and for Matlab by Alan Genz. Minor bug modifications 7/98, 2/10.
%
  if dh == inf || dk == inf, p = 0;
  elseif dh == -inf, if dk == -inf, p = 1; else p = Phi(-dk); end
  elseif dk == -inf, p = Phi(-dh);
  else
    if abs(r) < 0.3, ng = 1; lg = 3;
      %       Gauss Legendre points and weights, n =  6
      w(1:3,1) = [0.1713244923791705 0.3607615730481384 0.4679139345726904]';
      x(1:3,1) = [0.9324695142031522 0.6612093864662647 0.2386191860831970]';
    elseif abs(r) < 0.75,  ng = 2; lg = 6;
      %       Gauss Legendre points and weights, n = 12
      w(1:3,2) = [.04717533638651177 0.1069393259953183 0.1600783285433464]';
      w(4:6,2) = [0.2031674267230659 0.2334925365383547 0.2491470458134029]';
      x(1:3,2) = [0.9815606342467191 0.9041172563704750 0.7699026741943050]';
      x(4:6,2) = [0.5873179542866171 0.3678314989981802 0.1252334085114692]';
    else, ng = 3; lg = 10;
      %       Gauss Legendre points and weights, n = 20
      w(1:3,3) = [.01761400713915212 .04060142980038694 .06267204833410906]';
      w(4:6,3) = [.08327674157670475 0.1019301198172404 0.1181945319615184]';
      w(7:9,3) = [0.1316886384491766 0.1420961093183821 0.1491729864726037]';
      w(10,3) = 0.1527533871307259;
      x(1:3,3) = [0.9931285991850949 0.9639719272779138 0.9122344282513259]';
      x(4:6,3) = [0.8391169718222188 0.7463319064601508 0.6360536807265150]';
      x(7:9,3) = [0.5108670019508271 0.3737060887154196 0.2277858511416451]';
      x(10,3) = 0.07652652113349733;
    end, h = dh; k = dk; hk = h*k; bvn = 0;
    if abs(r) < 0.925, hs = ( h*h + k*k )/2; asr = asin(r);  
      for i = 1 : lg
	sn = sin( asr*( 1 - x(i,ng) )/2 );
	bvn = bvn + w(i,ng)*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
	sn = sin( asr*( 1 + x(i,ng) )/2 );
	bvn = bvn + w(i,ng)*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
      end, bvn = bvn*asr/( 4*pi );  
      bvn = bvn + Phi(-h)*Phi(-k);  
    else, twopi = 2*pi; if r < 0, k = -k; hk = -hk; end
      if abs(r) < 1, as = (1-r)*(1+r); a = sqrt(as); bs = (h-k)^2;
	c = ( 4 - hk )/8 ; d = ( 12 - hk )/16; asr = -( bs/as + hk )/2;
	if asr > -100  
	  bvn = a*exp(asr)*( 1 - c*(bs-as)*(1-d*bs/5)/3 + c*d*as*as/5 );
	end
	if hk > -100, b = sqrt(bs); sp = sqrt(twopi)*Phi(-b/a);
	  bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*( 1 - d*bs/5 )/3 );
	end, a = a/2;
	for i = 1 : lg
	  for is = -1 : 2 : 1, xs = ( a + a*is*x(i,ng) )^2; 
	    rs = sqrt( 1 - xs ); asr = -( bs/xs + hk )/2;
	    if asr > -100, sp = ( 1 + c*xs*( 1 + d*xs ) );
	      ep = exp( -hk*xs/( 2*(1+rs)^2 ) )/rs;
	      bvn = bvn + a*w(i,ng)*exp(asr)*( ep - sp );
	    end
	  end
	end, bvn = -bvn/twopi;
      end
      if r > 0, bvn =  bvn + Phi( -max( h, k ) ); 
      elseif h >= k, bvn = -bvn;
      else, if h < 0, L = Phi(k) - Phi(h); else, L = Phi(-h) - Phi(-k); end
        bvn =  L - bvn;
      end      
    end, p = max( 0, min( 1, bvn ) );
  end
%
%   end bvnu
%
