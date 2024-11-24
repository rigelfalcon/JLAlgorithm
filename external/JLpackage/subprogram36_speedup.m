


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

dtype=class(G);

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

Gm = zeros(K, K,dtype); % Preallocate

for i = 2:N
    for j = 1:N-i+1
        bInd = i+j-2; % Determine which block to use
        Gm((i-1)*M+1:i*M, (j-1)*M+1:j*M) =  G(:,M*bInd+1:M*(bInd+1));
    end
end

GmT = zeros(K, K,dtype); % Preallocate

for i = 2:N
    for j = 1:N-i+1
        bInd = i+j-2; % Determine which block to use
        GmT((i-1)*M+1:i*M, (j-1)*M+1:j*M) =  GT(:,M*bInd+1:M*(bInd+1));
    end
end
%disp(Gm);

Dt=zeros(2*K,dtype);
Dt(1:K,1:K)=D;   Dt(1:K,K+1:end)=-GmT;
Dt(K+1:end,1:K)=Gm;   Dt(K+1:end,K+1:end)=D;

%disp(Dt);    % det(Dt)
B1=zeros(2*K,M,dtype);   B1(1:M,1:M)=Im;
B2=zeros(2*K,M,dtype);   B2(K+1:K+M,1:M)=Im;
B=cat(2,B1,B2);

if M>128
    Dt=gpuArray(Dt);
    B=gpuArray(B);
end
X=Dt\B;
%X=Dt\B1;

if M>128
    X=gather(X);
end
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





%{
function U = subprogram36_speedup(G)
% Optimized version of subprogram36

[M, total_cols] = size(G);
N = total_cols / M;
K = M * N;

Im = eye(M);

% Compute GT by transposing blocks of G
GT = reshape(permute(reshape(G, M, M, N), [2, 1, 3]), M, M * N);

% Construct sparse D
diag_indices = 1:K;
row_indices = diag_indices;
col_indices = diag_indices;
values = ones(size(diag_indices));

% Additional indices for D
[kk, mm] = ndgrid(1:N-1, 1:M);
additional_row_indices = mm(:);
additional_col_indices = M * kk(:)' + mm(:)';
additional_values = ones(size(additional_row_indices));

% Combine indices
row_indices = [row_indices, additional_row_indices'];
col_indices = [col_indices, additional_col_indices'];
values = [values, additional_values'];

% Create sparse D
D = sparse(row_indices, col_indices, values, K, K);

% Reshape G and GT into blocks
G_blocks = reshape(G, M, M, N);
GT_blocks = reshape(GT, M, M, N);

% Generate i_vals, j_vals, bInd_vals
num_blocks = (N - 1) * N / 2;
i_vals = [];
j_vals = [];
bInd_vals = [];

for i = 2:N
    j_range = 1:N - i + 1;
    i_vals = [i_vals; repmat(i, length(j_range), 1)];
    j_vals = [j_vals; j_range'];
    bInd_vals = [bInd_vals; (i + j_range - 2)'];
end

% Prepare grid indices within a block
[mm, nn] = ndgrid(0:M-1, 0:M-1);
mm = mm(:); % size M^2 x 1
nn = nn(:); % size M^2 x 1

% Generate row_starts and col_starts
row_starts = (i_vals - 1) * M;
col_starts = (j_vals - 1) * M;

% Total elements
total_elements = M^2 * numel(i_vals);

% Generate row and column indices for Gm
Gm_row_indices = kron(row_starts, ones(M^2, 1)) + repmat(mm, numel(i_vals), 1) + 1;
Gm_col_indices = kron(col_starts, ones(M^2, 1)) + repmat(nn, numel(i_vals), 1) + 1;

% Extract blocks from G_blocks
block_indices = bInd_vals + 1;
Gm_blocks = G_blocks(:, :, block_indices);

% Flatten Gm_values
Gm_values = reshape(Gm_blocks, [], 1);

% Create sparse Gm
Gm = sparse(Gm_row_indices, Gm_col_indices, Gm_values, K, K);

% Do the same for GmT
GmT_blocks = GT_blocks(:, :, block_indices);
GmT_values = reshape(GmT_blocks, [], 1);
GmT = sparse(Gm_row_indices, Gm_col_indices, GmT_values, K, K);

% Construct Dt as sparse matrix
Dt_upper_left = D;
Dt_upper_right = -GmT;
Dt_lower_left = Gm;
Dt_lower_right = D;

Dt = [Dt_upper_left, Dt_upper_right; Dt_lower_left, Dt_lower_right];

% Construct B
B1 = sparse(2*K, M);
B1(1:M, 1:M) = Im;

B2 = sparse(2*K, M);
B2(K+1:K+M, 1:M) = Im;

B = [B1, B2];

% Solve the linear system
X = Dt \ B;

% Reshape X to get X11, X12, X21, X22
X11 = permute(reshape(X(1:K, 1:M), [M, N, M]), [1, 3, 2]);
X21 = permute(reshape(X(K+1:end, 1:M), [M, N, M]), [1, 3, 2]);
fX21 = flip(X21, 3);

X12 = permute(reshape(X(1:K, M+1:2*M), [M, N, M]), [1, 3, 2]);
X22 = permute(reshape(X(K+1:end, M+1:2*M), [M, N, M]), [1, 3, 2]);
fX22 = flip(X22, 3);

% Construct U
U = [cat(1, X11, fX21), cat(1, X12, fX22)];

end
%}