%01-05-2018: I wrote this program of finding the inverse of the polynomial
%when I find that inverse cannot be find by ifft when polynomial has zero
%on the boundary. The polynomial and its inverse have the same number of
%coefficients

function f_m=inverse_polynomial(b)

m=length(b)-1;

b0=b(1);

b=b./b0;

Q=1;
 B=[];

for n=1:m;
    B=[b(n+1) B];
   c=(-B*Q');
    Q=[Q c];     
end

f_m=Q/b0;

clear m