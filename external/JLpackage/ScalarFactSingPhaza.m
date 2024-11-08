% 01.05.2018: this program designed for singular matrices
function phasa=ScalarFactSingPhaza(a)   %R+1 entries are supplied

R=size(a,2)-1;
FTTP=2*R;

%e=FTTP^(-1);
e=10^(-4);

S=a;
%S=a.^2;
S=[S fliplr(S)];
S(R+2)=[];
S(end)=[];



tmp= e > a;
a(tmp)=2*e;

a=log(a);
%a=2*log(a);
aa=a;

a(1)=[];
a(end)=[];

aa=[aa fliplr(a)];

clear a

aa=ifft(aa);

aa=[aa(1) 2*aa(2:R) aa(R+1) zeros(1,R-1)];
%aa=[0.5*aa(1) aa(2:R) 0.5*aa(R+1) zeros(1,R-1)];

aa=fft(aa);
phasa=imag(aa);
% aa=i*imag(aa);
% phasa=exp(aa);




    