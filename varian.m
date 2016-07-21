function v = varian( H, b, h, k )
% this function computes the variances for 1,2,...,k-way random walk
I = eye(size(H,1));
v = zeros(k,1);
x = linsolve(I-H,b);
ex = x'*h;
h = h/abs(ex);
ex = 1;

hh = sum(abs(h))*abs(h);

[HH,~] = hhead(H,k);

TH = I;
for i=1:k
    TH = HH(:,:,i)*TH;
    G = I;
    if i>1
        temp = I;
        for j=1:i-1
            G = G+temp*HH(:,:,i-j+1);
            temp = temp*HH(:,:,i-j+1);
        end
    end

    if i==1
        if max(abs(eig(TH)))>=1
            fprintf('standard mc diverges\n');
            break
        end
    end
    bb = G*diag(b)*(2*H*x+b);
    t = linsolve(I - TH,bb);
    v(i) = hh'*(t) - ex*ex;
end