 function permsall = getComb(n,dim);
% getComb: all permutations. dim dimensions inside are ordered increasingly, 
% out dimensions are unordered. 
C = nchoosek(1:n,dim);

permsall = zeros(0,n);
for j=1:size(C,1)
    out = setdiff(1:n,C(j,:));
    po = perms(out);
    permsall = [permsall;ones(size(po,1),1)*C(j,:),po];
end;
