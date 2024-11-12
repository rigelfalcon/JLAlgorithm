function U = subprogram34_optimized(G)

[m, N] = size(G); B = G(m, :); I = eye(m);

ETA = zeros(N, m); ETA(1, m) = 1;

% Construct Dinv
c = zeros(1, 2*N-1);
c(1) = B(1);
r = [B zeros(1, N-1)];
Dinv = toeplitz(c(1:N), r(1:N));
Dinv(1, :) = -sum(Dinv(2:end, :), 1);
Dinv(1, 1) = 1;

% Preallocate Dt and Tsave
Dt = zeros(N);
Tsave = zeros(N, m-1, N);

% Precompute FFTs
Dinv_fft = fft([Dinv; zeros(N, N)], [], 1);
gm_fft_storage = zeros(2*N, m-1, N);

for k = 1:m-1
    Bk = flip(G(k, :))';
    c_k = zeros(N, 1);
    c_k(1) = Bk(end);
    gm_k = hankel(Bk, c_k);
    gm_k(1, :) = 0;

    % Compute P using FFTs
    gm_k_padded = [gm_k; zeros(N, N)];
    gm_fft = fft(gm_k_padded, [], 1);
    gm_fft_storage(:, k, :) = gm_fft;

    P_fft = Dinv_fft .* gm_fft;
    P = ifft(P_fft, [], 1);
    P = P(1:N, :);

    % Update Tsave and ETA
    Tsave(:, k, :) = P;
    ETA(:, k) = -P(:, 1);

    % Update Dt using efficient computation
    Dt = Dt + P * P;
end

Dt = Dt + eye(N);
X = Dt \ ETA;

% Construct U
U = zeros(m, m, N);
U(m, :, :) = permute(X, [2, 3, 1]);

% Compute Uabove using Tsave and X
Uabove = zeros(m-1, m, N);
for k = 1:m-1
    P = squeeze(Tsave(:, k, :));
    Uabove(k, :, :) = P * X;
end

U(1:m-1, :, :) = Uabove;
I(m, m) = 0; U(:, :, 1) = U(:, :, 1) + I;

end
