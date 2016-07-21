function [HH,lam] = hhead(H,k)
% this function compute the tensor HH (\head{H})
% and the coresponding lambda values
lam = zeros(k,1);
n = size(H,1);
e = ones(n,1);
absh = abs(H);

w = e;
HH = zeros(n,n,k);
for i=1:k
    eta = absh*w;
    HH(:,:,i) = diag(eta)*absh*diag(1./w);
    w = eta;
    lam(i) = (max(eta))^2;
end