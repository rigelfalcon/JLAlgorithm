
function f_inv_phi=inverse_polynomial_block(phi,f)

[M,~,N]=size(phi);

f_inv_phi=zeros(M,M,N);


f_inv_phi(:,:,1)=f(:,:,1)\phi(:,:,1);
for k=2:N
    S=zeros(M);
    for j=1:k-1
        S=S+f(:,:,k-j+1)*f_inv_phi(:,:,j);
        f_inv_phi(:,:,k)=f(:,:,1)\(phi(:,:,k)-S);
    end
end

end