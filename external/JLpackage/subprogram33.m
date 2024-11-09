% This subprogram has been changed in order to fit this JLE-4_Cholesky

function U=subprogram33(G)
% N, m positive integers, m>=2
% B is a column vector N+1, G is a matrix Nx(m-1)
G=G.';
[N,m]=size(G);

B=G(:,m);
I=eye(m,class(G));

%%

%[X,ETA]=subprogram2(B,G);
% %{
ETA=zeros(N,m,class(G));

for i=1:m-1
    eta=conv(B,G(:,i));
    eta=flipud(eta(1:(N)));
    ETA(:,i)=eta;
end;
% %}

%%
%{
ETA = zeros(N, m, class(G));
eta= convn(G(:,1:m-1),B);
ETA(:,1:m-1)=eta(1:N, :);
%}
%%
%{
% Preallocate and compute ETA using vectorized FFT-based convolution
NN = 2*N - 1;
B_fft = fft(B, NN);
G_fft = fft(G(:,1:m-1), NN);
ETA_fft = B_fft .* G_fft;
eta = ifft(ETA_fft, NN);
eta = flipud(eta(1:N, :));
ETA = zeros(N, m, class(G));
ETA(:,1:m-1) = eta;
% error=sum((ETA2-ETA).^2,'all')
%}
%%
% Initialize G0 and other variables
G0=ETA;

clear B  G

ETA(1,m)=1;  %this is a new here

G0(N,m)=1;

U=zeros(N,N,class(ETA));

D=zeros(1,N,class(ETA));

% for j=1:(N-1)
%     g0=G0(N+1-j,:);
%     U0=G0*g0';
%     d0=U0(end);
%     D(N-j+1)=d0;
%     % V0=U0(2:(end))-U0(1:(end-1));
%     V0 = diff(U0);
%     M0=V0*g0/d0;
%     G0=M0+G0(1:(end-1),:);
%     U(1:(N+1-j),N+1-j)=U0;
% 
% 
% end;

G0_size = N;
for j = 1:(N-1)
    idx = N+1-j;
    g0 = G0(idx, :);
    U0 = G0(1:G0_size, :) * g0';
    d0 = U0(end);
    D(N-j+1) = d0;
    V0 = diff(U0);
    M0 = (V0 * g0) / d0;
    G0 = M0 + G0(1:G0_size-1, :);
    U(1:G0_size, idx) = U0;
    G0_size = G0_size - 1;
end;

g0=G0(1,:);
U0=G0*g0';
d0=U0(end);
D(1)=d0;
U(1,1)=U0;



clear G0

Y=U\ETA;




Y=diag(D)*Y;
X=U.'\Y;

ETA(:,m)=[];

clear U  U1



U=zeros(m,m,N,class(ETA));

U(m,:,:)=X.';

X=flip(X,1);  %this has been changed from Subprogram 3 (moved down); U keeps the last row without "flip"




NN=2*N-1;  ETA=ETA.';
ETA=fft(ETA,NN,2);  X=fft(X,NN,1);
% pp=zeros(m-1,m,NN);


% for i=1:NN
%     pp(:,:,i)=ETA(:,i)*X(i,:);
% end
ETA=permute(ETA,[1,3,2]);
X=permute(X,[3,2,1]);
pp=ETA.*X;

pp=ifft(pp,NN,3);

U(1:m-1,:,1:N)=pp(:,:, N:NN);


% for j=1:m
%     for i=1:(m-1)
%         pp=conv(X(:,j),ETA(:,i));
%         U(i,j,:)=pp((N+1):(2*N+1));
%     end;
% end;

I(m,m)=0;   U(:,:,1)=U(:,:,1)-I;

v=sum(U,3);
%u=inv(v);

% for k=1:N
%     U(:,:,k)=U(:,:,k)/v;
% end;
U=pagemrdivide(U,v);

%----------here we check the unitarity of U (flipping the last row)
%         UT=U;
%         for k=1:m
%            UT(m,k,:)=flip(squeeze(UT(m,k,:)));
%         end
%         UT_flip=conj(permute(UT,[2,1,3]));
%   UT_flip=flip(UT_flip,3);
%
%   Id=Prod_Mat_Pol(UT,UT_flip);
%      II=zeros(m,m,2*N-1);
%      II(:,:,N)=eye(m);
%      ErrUnit=max(max(max(abs(Id-II))))
% %
end