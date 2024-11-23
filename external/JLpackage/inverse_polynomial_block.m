function f_inv_phi = inverse_polynomial_block(phi, f)

[M, ~, N] = size(phi);

f_inv_phi = zeros(M, M, N);

% Precompute the inverse of f(:,:,1)
inv_f1 = inv(f(:, :, 1));

% Compute the first term
f_inv_phi(:, :, 1) = inv_f1 * phi(:, :, 1);

for k = 2:N
    % Indices for the current convolution
    idx_f = k:-1:2;          % Indices for f matrices
    idx_phi = 1:k-1;         % Indices for f_inv_phi matrices

    % Extract the required slices
    f_matrices = f(:, :, idx_f);                % [M, M, k - 1]
    f_inv_phi_matrices = f_inv_phi(:, :, idx_phi); % [M, M, k - 1]

    % Perform batch matrix multiplication
    products = pagemtimes(f_matrices, f_inv_phi_matrices); % [M, M, k - 1]

    % Sum over the third dimension to compute S
    S = sum(products, 3);

    % Compute f_inv_phi(:, :, k)
    f_inv_phi(:, :, k) = inv_f1 * (phi(:, :, k) - S);
end

end





% function f_inv_phi=inverse_polynomial_block(phi,f)
% 
% [M,~,N]=size(phi);
% 
% f_inv_phi=zeros(M,M,N);
% 
% 
% f_inv_phi(:,:,1)=f(:,:,1)\phi(:,:,1);
% for k=2:N
%     S=zeros(M);
%     for j=1:k-1
%         S=S+f(:,:,k-j+1)*f_inv_phi(:,:,j);
%         f_inv_phi(:,:,k)=f(:,:,1)\(phi(:,:,k)-S);
%     end
% end
% 
% end