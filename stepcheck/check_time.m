% r=1000;%00
% n=9;
% NFFT=2^nextpow2(n);
% U=randn(1,1,NFFT);
% tic
% for ir=2:r
%     U=padarray(U,[1,1],0,'post');
%     U(end,end,:)=1;
%     Unew=randn(ir,ir,NFFT);
% end
% toc

r=1000;%00
n=9;
NFFT=2^nextpow2(n);
U=zeros(r,r,NFFT);
U(1,1,NFFT)=randn(1,1,NFFT);
tic
for ir=2:r
    % U=padarray(U,[1,1],0,'post');
    U(end,end,:)=1;
    Unew=randn(ir,ir,NFFT);
end
toc

% r=1000;%00
% n=9;
% NFFT=2^nextpow2(n);
% U=randn(1,1,NFFT);
% tic
% for ir=2:r
%     U=padarray(U,[1,1],0,'post');
%     U(end,end,:)=1;
%     Unew=randn(ir,ir,NFFT);
% end
% toc