clear
M=5;
N=10;

% b=randn(M,M*N);
% f_m=inverse_polynomial_matrix(b);

phi=randn(M,M,N);
f=randn(M,M,N);
% f_inv_phi=inverse_polynomial_block(phi,f);

f_inv_phi=zeros(M,M,N);


f_inv_phi(:,:,1)=f(:,:,1)\phi(:,:,1);
for k=2:N
    S=zeros(M);
    for j=1:k-1
        S=S+f(:,:,k-j+1)*f_inv_phi(:,:,j);
        f_inv_phi(:,:,k)=f(:,:,1)\(phi(:,:,k)-S);
    end
end


phi2=Prod_Mat_Pol(f,f_inv_phi);

phi2=phi2(:,:,1:N);

diff=abs(phi2-phi);
max(abs(phi2-phi),[],'all')





