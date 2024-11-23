

function U=subprogram36_speedup(G)
% clc

% N=3;   M=2;
%N=10;  M=7;
%
% g1=[1 2 3;  4 5 6; 7 8 9];
% g2=[10 11 12; 13 14 15; 16 17 18];
% g3=[19 20 21; 22 23 24; 25 26 27];
% g4=[28 29 30; 31 32 33; 34 35 36];
[M,N]=size(G);      N=N/M;      K=M*N;



Im=eye(M);

%G=zeros(M*2,N*M);
% G=zeros(M,N*M);
% G=rand(M,N*M);

% G(:,M+1:2*M)=g1;   G(:,2*M+1:3*M)=g2;
% G(:,3*M+1:4*M)=g3;   G(:,4*M+1:5*M)=g4;

GT=G;
% GT=reshape(GT, K, K, M);
% GT=flip(GT,3);
% GT=reshape(GT, K, K*M)
GT = reshape(permute(reshape(GT, M, M, N), [2, 1, 3]), M, M * N);






D = eye(K, K); % Preallocate
for k=1:N-1
    D(1:M,M*k+1:M*(k+1))=Im;
end

% D;

Gm = zeros(K, K); % Preallocate

for i = 2:N
    for j = 1:N-i+1
        bInd = i+j-2; % Determine which block to use
        Gm((i-1)*M+1:i*M, (j-1)*M+1:j*M) =  G(:,M*bInd+1:M*(bInd+1));
    end
end


GmT = zeros(K, K); % Preallocate

for i = 2:N
    for j = 1:N-i+1
        bInd = i+j-2; % Determine which block to use
        GmT((i-1)*M+1:i*M, (j-1)*M+1:j*M) =  GT(:,M*bInd+1:M*(bInd+1));
    end
end

% disp(Gm);
% Gm=fast_hankel(G');
% Gm(1,:)=0;
% GT=fast_hankel(GT');
% Gm(1,:)=0;


Dt=zeros(2*K);
Dt(1:K,1:K)=D;   Dt(1:K,K+1:end)=-GmT;
Dt(K+1:end,1:K)=Gm;   Dt(K+1:end,K+1:end)=D;

%disp(Dt);    % det(Dt)
B1=zeros(2*K,M);   B1(1:M,1:M)=Im;
B2=zeros(2*K,M);   B2(K+1:K+M,1:M)=Im;
B=cat(2,B1,B2);

X=Dt\B;
%X=Dt\B1;

X11=permute(reshape(X(1:K,1:M),[M,N,M]),[1,3,2]);

X21=permute(reshape(X(K+1:end,1:M),[M,N,M]),[1,3,2]);
fX21=flip(X21,3);

X12=permute(reshape(X(1:K,M+1:2*M),[M,N,M]),[1,3,2]);

X22=permute(reshape(X(K+1:end,M+1:2*M),[M,N,M]),[1,3,2]);
fX22=flip(X22,3);


U=cat(2,cat(1,X11,fX21),...
    cat(1,X12,fX22));

% UT=U;
% 
% size(U);   %  6x6x5
% 
% for j=1:size(U,3)
%     UT(:,:,j)=UT(:,:,j)';
% end
% UT=flip(UT,3);
% PrM=Prod_Mat_Pol(U,UT)


end


function [Gm, GmT] = optimize_Gm_GmT_intermediate(G, GT, M, N)

% Calculate K
K = M * (N - 1);

% Preallocate Gm and GmT
Gm = zeros(K, K);
GmT = zeros(K, K);

for i = 2:N
    j_range = 1:N - i + 1;
    bInd_range = i + j_range - 2;
    
    row_idx = (i - 1) * M + 1:i * M;
    col_starts = (j_range - 1) * M + 1;
    G_col_starts = bInd_range * M + 1;
    
    % Vectorize over j_range
    for idx = 1:length(j_range)
        col_idx = col_starts(idx):col_starts(idx) + M - 1;
        G_col_idx = G_col_starts(idx):G_col_starts(idx) + M - 1;
        
        % Assign blocks
        Gm(row_idx, col_idx) = G(:, G_col_idx);
        GmT(row_idx, col_idx) = GT(:, G_col_idx);
    end
end

end


% function [Gm, GmT] = optimize_Gm_GmT(G, GT, M, N)
% 
% % Calculate K
% K = M * (N - 1);
% 
% % Preallocate Gm and GmT
% Gm = zeros(K, K);
% GmT = zeros(K, K);
% 
% % Precompute total number of blocks
% total_blocks = N * (N - 1) / 2;
% 
% % Initialize arrays to store indices
% row_indices = zeros(M, total_blocks);
% col_indices = zeros(M, total_blocks);
% G_indices = zeros(M, total_blocks);
% G_indices_T = zeros(M, total_blocks);
% 
% % Initialize block counter
% block_counter = 1;
% 
% % Generate indices
% for i = 2:N
%     j_range = 1:N - i + 1;
%     bInd_range = i + j_range - 2;
% 
%     row_start = (i - 1) * M + 1;
%     row_end = i * M;
%     row_idx = row_start:row_end;
% 
%     for idx = 1:length(j_range)
%         j = j_range(idx);
%         bInd = bInd_range(idx);
% 
%         col_start = (j - 1) * M + 1;
%         col_end = j * M;
%         col_idx = col_start:col_end;
% 
%         G_col_start = bInd * M + 1;
%         G_col_end = (bInd + 1) * M;
%         G_col_idx = G_col_start:G_col_end;
% 
%         % Store indices
%         row_indices(:, block_counter) = row_idx;
%         col_indices(:, block_counter) = col_idx;
%         G_indices(:, block_counter) = G_col_idx;
%         G_indices_T(:, block_counter) = G_col_idx; % Assuming GT has the same indexing
% 
%         block_counter = block_counter + 1;
%     end
% end
% 
% % Flatten the indices
% row_indices_flat = row_indices(:);
% col_indices_flat = col_indices(:);
% G_indices_flat = G_indices(:);
% G_indices_T_flat = G_indices_T(:);
% 
% % Use linear indexing to assign values
% Gm(sub2ind(size(Gm), repmat(row_indices_flat, M, 1), reshape(repmat(col_indices_flat', M, 1), [], 1))) = ...
%     G(repmat(1:M, total_blocks * M, 1), G_indices_flat);
% 
% GmT(sub2ind(size(GmT), repmat(row_indices_flat, M, 1), reshape(repmat(col_indices_flat', M, 1), [], 1))) = ...
%     GT(repmat(1:M, total_blocks * M, 1), G_indices_T_flat);
% 
% end
