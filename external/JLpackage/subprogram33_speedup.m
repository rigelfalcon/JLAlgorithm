function U = subprogram33_speedup(G)
    G = G.';
    [N, m] = size(G);
    B = G(:, m);
    
    % Preallocate and compute ETA using vectorized FFT-based convolution
    NN = 2*N - 1;
    B_fft = fft(B, NN);
    G_fft = fft(G(:,1:m-1), NN);
    ETA_fft = B_fft .* G_fft;
    eta = ifft(ETA_fft, NN);
    eta = flipud(eta(1:N, :));
    ETA = zeros(N, m, class(G));
    ETA(:,1:m-1) = eta;
    
    % Initialize G0 and other variables
    ETA(1, m) = 1;
    G0 = ETA;
    G0(N, m) = 1;
    U = zeros(N, N, class(ETA));
    D = zeros(1, N, class(ETA));
    
    % Optimize the loop updating G0 and U
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
    g0 = G0(1, :);
    U0 = G0 * g0';
    d0 = U0(end);
    D(1) = d0;
    U(1,1) = U0;
    
    % Solve linear systems using U
    Y = U \ ETA;
    Y = bsxfun(@times, D.', Y);
    X = U.' \ Y;
    ETA(:, m) = [];
    
    % % Build U matrix using optimized FFT operations
    % U = zeros(m, m, N, class(ETA));
    % U(m,:,:) = permute(X, [3, 2, 1]);
    % X = flip(X, 1);
    % ETA_fft = fft(ETA.', NN, 2);
    % X_fft = fft(X, NN, 1);
    % pp=ifft(ETA_fft.*X_fft,NN,3);
    % % pp = ifft(bsxfun(@times, ETA_fft, permute(X_fft, [3,2,1])), NN, 3);
    % U(1:m-1,:,1:N) = pp(:,:, N:NN);
    U=zeros(m,m,N,class(ETA));

    U(m,:,:)=X.';
    
    X=flip(X,1);  %this has been changed from Subprogram 3 (moved down); U keeps the last row without "flip"

    NN=2*N-1;  ETA=ETA.';
    ETA=fft(ETA,NN,2);  X=fft(X,NN,1);

    ETA=permute(ETA,[1,3,2]);
    X=permute(X,[3,2,1]);
    pp=ETA.*X;
    pp=ifft(pp,NN,3);
    U(1:m-1,:,1:N)=pp(:,:, N:NN);


    % Final adjustments and normalization
    I = eye(m, class(U)); I(m,m) = 0;
    U(:,:,1) = U(:,:,1) - I;
    v = sum(U, 3);
    U = pagemrdivide(U, v);
end
