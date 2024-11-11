
d=80;
n=9;
FFTP=2^n;
nepoch=10;
% randn(d,d,FFTP)
x=randn(d,FFTP,nepoch);

X=fft(x,[],2);
X=permute(X,[1,4,2,3]);
% S=pagemtimes(X,'none',X,'ctranspose');
% S=
S=mean(S,4);






