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

%AA=rand(d,d,n+1);
 %AA=-1+2*AA;    %load('AA.mat');  %takes values in [-1;1]
 AA_flip=permute(AA,[2,1,3]);
 AA_flip=flip(AA_flip,3);
 A=Prod_Mat_Pol(AA,AA_flip); 
 
 tic
I=eye(d);
n_mid=(size(A,3)+1)/2;      %this part copied from Test_Wilson_Iter
S_ext=cat(3, A(:,:,n_mid:end), zeros(d,d, FFTP-2*n_mid+1), A(:,:,1:n_mid-1));
   S_ext=fft(S_ext,[],3);      %cheked positive definite
   S_ext=(S_ext+conj(permute(S_ext,[2 1 3])))/2;
%    for k=1:FFTP                         %We add ep on the diagonal to Wilson
%        S_ext(:,:,k)=S_ext(:,:,k);%+epp*I;
%    end
%    
   
   %--------------------------cpoied from GenerationMukesh
   tic
   
 A_ext=zeros(d,d,FFTP);
 
 A_ext(:,:,1)=chol(S_ext(:,:,1), 'lower');
 
   
    A_ext(:,:,FF2+1)=chol( S_ext(:,:,FF2+1), 'lower');
    
    for k=2:FF2
    A_ext(:,:,k)=chol(S_ext(:,:,k), 'lower');
    A_ext(:,:,FFTP+2-k)=conj(A_ext(:,:,k));
    end

Time_Cholesky=toc;

 %------------------here ends Cholesky factorization 
  
 Phaza=zeros(d,FFTP); 
 for r=1:d
 Diag=ScalarFactSingPhaza(squeeze(A_ext(r,r,1:FF2+1)).');
 Phaza(r,:)=exp(1j*Diag);
 end
 
 clear Diag
 
 

for k=1:FFTP
A_ext(:,:,k)=A_ext(:,:,k)*diag(Phaza(:,k));
end

clear Phaza

%--------here ends making analitic diagonal and changing phazes of off-diagonal terms  



for r=2:d
  r
  
  phi_temp=reshape(A_ext(r,1:r,:), [r,FFTP]);
    
%      phi_temp(r,:)=1./phi_temp(r,:);
%     
    phi_temp=ifft(phi_temp,FFTP,2);
    
     if r==2
      phi_max=abs(phi_temp(1:r-1,FF2:end));
      else
      phi_max=max(abs(phi_temp(1:r-1,FF2:end)));
       end
       temp=phi_max>10^(-4);
       
N = FF2-find(temp~=0, 1, 'first');
    
    b=phi_temp(r,1:N+1);    
     b=inverse_polynomial(b);
    
    G=zeros(r,N+1);
    
      G(r,:)=b(1,:);
     
      
        G(1:r-1,1:N)=phi_temp(1:r-1,FFTP-N+1:FFTP);
   
       
    clear  phi_temp
    
    U=subprogram33(G);
  
     U_ext=fft(U,FFTP,3);
    
    U_ext(r,:,:)=conj(U_ext(r,:,:));
    
    for k=1:FFTP
      A_ext(:,1:r,k)=A_ext(:,1:r,k)*U_ext(:,:,k);
    end
       
end
toc
%--------------HERE WE CHECK THE FINAL RESULT
%---- we check error in the frequency domain

   A=zeros(d,d,FF2);
    for k=1:FF2
    A(:,:,k)=A_ext(:,:,k)*A_ext(:,:,k)';
    end
    
    B=A-S_ext(:,:,1:FF2);     
 Error_frequency_domain=max(max(max(abs(B))))
 
 %--- we check if the result is causal (half of Fourier coefficients in the time domain should be close to 0)
 
 A_time_domain=ifft(A_ext,FFTP,3);
 ans1=A_time_domain(:,:,FF2+1:FFTP);           
  Check_causality=max(max(max(abs(ans1))))


 %-- Here we write the result in the form where S(0) is positive definite

cov1=A_time_domain(:,:,1);
Cov=cov1*cov1';   Cov=sqrtm(Cov);
cov2=inv(cov1); cov2=cov2*Cov;
%Det=det(Cov) 
for k=1:n+1
    A_time_domain(:,:,k)=A_time_domain(:,:,k)*cov2;
end

  Pasuxi_22=A_time_domain(:,:,1:n+1);