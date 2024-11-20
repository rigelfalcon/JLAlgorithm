%{
clc

N=3;   M=2;
%N=10;  M=7;

g1=[1 2; 3 4]; g2=[5 6 ; 7 15 ];  
d0=[9 10; 11 12]; d1=[13 14; 15 16]; d2=[17 18; 19 20];

Im=eye(M);

G=zeros(M*2,N*M);
% G=rand(M*2,N*M);

G(1:M,M+1:2*M)=g1;   G(1:M,2*M+1:3*M)=g2;
G(M+1:2*M,1:M)=d0; G(M+1:2*M,M+1:2*M)=d1;   G(M+1:2*M,2*M+1:3*M)=d2;  


G1=G(1:M,:);
G2=G(M+1:2*M,:);
%}


function U=subprogram35(G)
[m,N]=size(G);   M=m/2;   N=N/M;      K=M*N;
G1=G(1:M,:);
G2=G(M+1:2*M,:);
Im=eye(M);

D = zeros(K, K); % Preallocate

for i = 2:N
    for j = i:N
        bInd = abs(i-j); % Determine which block to use
        D((i-1)*M+1:i*M, (j-1)*M+1:j*M) =  G2(:,M*bInd+1:M*(bInd+1));
    end
end

D1=repmat(Im, 1, N);
D(1:M,:)=D1;
% disp(D);

Gm = zeros(K, K); % Preallocate

for i = 2:N
    for j = 1:N-i+1
        bInd = i+j-2; % Determine which block to use
        Gm((i-1)*M+1:i*M, (j-1)*M+1:j*M) =  G1(:,M*bInd+1:M*(bInd+1));
    end
end
% disp(Gm);


Dt=zeros(2*K);
Dt(1:K,1:K)=D;   Dt(1:K,K+1:end)=-Gm;
Dt(K+1:end,1:K)=Gm;   Dt(K+1:end,K+1:end)=D;

 % disp(Dt);    % det(Dt)
 B1=zeros(2*K,M);   B1(1:M,1:M)=Im;
 %B2=zeros(2*K,M);   B2(K+1:K+M,1:M)=Im;
 %B=cat(2,B1,B2);

 % X=Dt\B;
 X=Dt\B1;

 X1=permute(reshape(X(1:K,:),[M,N,M]),[1,3,2]);
 X2=permute(reshape(X(K+1:end,:),[M,N,M]),[1,3,2]);

 U=cat(1,cat(2,X1,-X2),...
         cat(2,X2,X1));

end



% X =
% 
%     0.8665    0.9293   -0.0207   -0.0752
%    -0.0985   -0.1039    0.0138    0.0500
%     0.5691    0.6459   -0.4711   -0.5185
%    -0.2956   -0.3214    0.0547    0.0820
%    -0.4356   -1.5752    0.4918    0.5937
%     0.3941    1.4253   -0.0685   -0.1320
%     0.0207    0.0752    0.8665    0.9293
%    -0.0138   -0.0500   -0.0985   -0.1039
%     0.4711    0.5185    0.5691    0.6459
%    -0.0547   -0.0820   -0.2956   -0.3214
%    -0.4918   -0.5937   -0.4356   -1.5752
%     0.0685    0.1320    0.3941    1.4253


