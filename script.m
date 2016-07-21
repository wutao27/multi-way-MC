%% -- synthetic experiments: compute success ratio for given RHO
RHO = 0.9;
count = zeros(5,1);
for i=1:100
    n = 1000;
    A = sprand(n,n,0.2);
    rho = max(abs(eig(full(abs(A)))));
    A = (RHO/rho)*A;
    [~,lam] = hhead(A,5);
    count = count + double(lam<1);
end
ratio = count/100;

%% - synthetic experiments: plot the success ratio
clear all
load('ratio')
for i=1:5
    plot(x,y(i,:),'-o','LineWidth',4)
    hold on
end
set(gca,'fontsize',20)
h_legend = legend('standard 1-way', '2-way', '3-way','4-way', '5-way');
set(h_legend,'FontSize',20);
xlabel('$\rho(\tilde{\mathbf{H}}^{+})$','Interpreter','latex','FontSize',24,'FontWeight','bold');
ylabel('ratio','FontSize',24,'FontWeight','bold')

%% - synthetic experiments: compute the speed-up ratio for given RHO
clear all
RHO = 0.8;
n = 1000;
result = zeros(5,1);
count = 1;
while count <= 100
    A = sprand(n,n,0.2);
    rho = max(abs(eig(full(abs(A)))));
    A = (RHO/rho)*A;
    b = rand(n,1);
    h = rand(n,1);
    v = varian( A, b, h, 5 );
    if v(1)>0
        result = result + v(1)./v;
        count = count+1;
    end
end
result = result/100;


%% - Harwell-Boeing sparse matrix collection
% compute the average speed-up ratio
clear all
load('hbset.mat')
result = zeros(5,1);
count = 0;
for i=1:size(hbset,1)
   mat=hbset(i);
   Prob = UFget(mat);
   A = Prob.A;
   n = size(A,1);
   d = diag(A);
   d(abs(d) < eps(1)) = 1;
   H = speye(n)-diag(1./d)*A;
%    deal with the cases when there are empty rows in H
   zeroind = find((abs(H)*ones(n,1))==0);
   while size(zeroind,1)~=0
    H(zeroind,:) = [];
    H(:,zeroind) = [];
    n = size(H,1);
    if n==0
       continue
    end
    zeroind = find((abs(H)*ones(n,1))==0);
   end
        
    b = rand(n,1);
    h = rand(n,1);
    v = varian( H, b, h, 5 );
     if v(1)>0 && abs(v(1)/v(2) - 1) > 0.001
         result = result + v(1)./v;
         count = count + 1;
     end
end
result = result/count;
result = result';
save hb_speed_ratio result
