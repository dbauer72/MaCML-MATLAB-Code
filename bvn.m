function p = bvn( xl, xu, yl, yu, r )
%BVN
%  A function for computing bivariate normal probabilities.
%  bvn calculates the probability that 
%    xl < x < xu and yl < y < yu, 
%  with correlation coefficient r.
%   p = bvn( xl, xu, yl, yu, r )
%

%   Author
%       Alan Genz, Department of Mathematics
%       Washington State University, Pullman, Wa 99164-3113
%       Email : alangenz@wsu.edu
%
  p = bvnu(xl,yl,r) - bvnu(xu,yl,r) - bvnu(xl,yu,r) + bvnu(xu,yu,r); 
  p = max( 0, min( p, 1 ) );
%
%   end bvn
%
