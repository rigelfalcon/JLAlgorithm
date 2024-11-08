function L = pagechol(A)
%PAGECHOL Page-wise Cholesky decomposition
%   L = PAGECHOL(A) performs the Cholesky decomposition on each page of A.
%   A can be an N-D array, where the first two dimensions correspond to
%   square matrices, and the remaining dimensions are pages.
%   The output L is an array of the same size as A, containing the lower
%   triangular Cholesky factors for each page.

% Check that the input array has at least 2 dimensions
if ndims(A) < 2
    error('Input A must be at least a 2-D array.');
end

% Get the size of the input array
sz = size(A);
nRows = sz(1);
nCols = sz(2);

% Check that the matrices are square
if nRows ~= nCols
    error('Matrices must be square in the first two dimensions.');
end

% Calculate the total number of matrices/pages
nMats = prod(sz(3:end));

% Reshape A to a 3-D array for easier processing
A = reshape(A, nRows, nCols, nMats);

% Preallocate the output array L
L = zeros(size(A));

% Loop over each matrix/page
for k = 1:nMats
    Ak = A(:, :, k);
    
    % Check if the matrix is symmetric
    if ~isequal(Ak, Ak')
        error('Matrix number %d is not symmetric.', k);
    end
    
    % Perform the Cholesky decomposition
    [Lk, p] = chol(Ak, 'lower');
    
    % Check if the matrix is positive definite
    if p ~= 0
        error('Matrix number %d is not positive definite.', k);
    end
    
    % Store the result in the output array
    L(:, :, k) = Lk;
end

% Reshape L back to the original size of A
L = reshape(L, sz);
end
