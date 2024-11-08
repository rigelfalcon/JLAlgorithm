clc
clear all
format long
d=100;        % matrix dimension
n=10;       % random polynomial degree
p=9;     FFTP=2^p;  % Number of FFT nodes in frequency domain
FF2=FFTP/2; % N=1000;

AA=zeros(d,d,n+1);

for k=1:n+1
    AD=randn(d,d); AA(:,:,k)=AD;
end


AA_flip=permute(AA,[2,1,3]);
AA_flip=flip(AA_flip,3);
A=Prod_Mat_Pol(AA,AA_flip);

tic
I=eye(d);
n_mid=(size(A,3)+1)/2;      %this part copied from Test_Wilson_Iter
S_ext=cat(3, A(:,:,n_mid:end), zeros(d,d, FFTP-2*n_mid+1), A(:,:,1:n_mid-1));
S_ext=fft(S_ext,[],3);      %cheked positive definite
S_ext=(S_ext+conj(permute(S_ext,[2 1 3])))/2;




