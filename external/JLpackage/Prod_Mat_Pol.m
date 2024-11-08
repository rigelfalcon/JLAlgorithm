function PrM=Prod_Mat_Pol(MP1,MP2)   

[d1,d2,m]=size(MP1);
m=m-1;
[d2,d3,n]=size(MP2);
n=n-1;
mn=m+n;

PrM=zeros(d1,d3,m+n+1);



for k=1:m+1
     temp=zeros(d1,d3,n+1); 
    for l=1:n+1
     temp(:,:,l)=MP1(:,:,k)*MP2(:,:,l); 
    end
PrM(:,:,k:k+n)=PrM(:,:,k:k+n)+temp;
end
end