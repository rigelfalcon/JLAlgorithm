% in this program we construct U without application of the displacement
% structure, solving directly the system Dt*X=B, where the first row of
% Dt is modified


function U=subprogram34(G)

[m,N]=size(G); B=G(m,:);   I=eye(m);

ETA=zeros(N,m);    ETA(1,m)=1;

% for i=1:m-1
%     eta=conv(B,G(:,i));
%     eta=flipud(eta(1:(N)));
%     ETA(:,i)=eta;
% end



c=zeros(1,N);

c(1)=B(1);   
Dinv=toeplitz(c,B);  

Dinv(1,:)=zeros(1,N);   s=sum(Dinv,1);  Dinv(1,:)=-s;  Dinv(1,1)=1;

Dt=zeros(N);

Tsave=zeros((m-1)*N,N);


f=0:N-1;
[fx,fy]=ndgrid(f);

%
B=flip(G(1,:))'; 
Bf=griddedInterpolant(f,B,'nearest','none');



for k=1:m-1
    B=flip(G(k,:))'; 
    % c(1)=B(end);
    % gm=hankel(B,c);
    % gm2=hankel(B);
    % Bf=griddedInterpolant(f,B,'nearest','none');
    Bf.Values=B;
    gm=Bf(fx+fy);

    gm(isnan(gm))=0;
    gm(1,:)=zeros;
    P=Dinv*gm;
    Tsave(N*(k-1)+1:N*(k-1)+N,:)=P;
    ETA(:,k)=-P(:,1);
    Dt=Dt+P*P;
end
%

% B=G(1:m-1,:)';
% Bf=griddedInterpolant(f,B,'nearest','none');
% gm=Bf(fx+fy);
% gm(isnan(gm))=0;
% gm=reshape(gm,[N,N,m-1]);
% P=pagemtimes(Dinv,gm);
% Tsave(1+1:N*(k-1)+N)



Dt=Dt+eye(N);
X=Dt\ETA;

% ETA(:,m)=[];


U=zeros(m,m,N);

U(m,:,:)=X.';  

Uabove=Tsave*X;

% Uabove=Uabove';
Uabove=reshape(Uabove,[N,m-1,m]);
Uabove=permute(Uabove,[2,3,1]);

% size(Uabove)





% % % % % X=flip(X,1);  %this has been changed from Subprogram 3 (moved down); U keeps the last row without "flip"
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % NN=2*N-1;  ETA=ETA.';   
% % % % % ETA=fft(ETA,NN,2);  X=fft(X,NN,1);
% % % % % pp=zeros(m-1,m,NN);
% % % % % 
% % % % % 
% % % % % for i=1:NN
% % % % %     pp(:,:,i)=ETA(:,i)*X(i,:);
% % % % % end
% % % % % 
% % % % % pp=ifft(pp,NN,3);
% % % % % 

% for k1=1:m-1
%     for k2=1:m
% % U(1:m-1,:,1:N)=Uabove;
% U(k1,k2,1:N)=squeeze(Uabove(N*(k1-1)+1:N*(k1-1)+N,k2));
%     end
% end

U(1:m-1,:,1:N)=Uabove;


I(m,m)=0;   U(:,:,1)=U(:,:,1)+I;