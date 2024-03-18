rng("default");
load("west0479");
load("datasets.mat");

A = dw8192;
m = 100;
tau = 1; % shift-and-invert spectral transformation

[n, ~] = size(A);
V = zeros(n, m+1);
H = zeros(m+1,m);

% Initialize first vector
V(:, 1) = rand([n, 1]);

% Normalize first vector
V(:, 1) = V(:, 1) / norm(V(:, 1));

% LU Factorization
[L, R, P] = lu(A - tau*speye(n));

% Get inverse  of A using LU factorization
inv_A = zeros(n, n);
B = eye(n);
for i = 1:n
    y = L \ (P * B(:, i)); % Forward substitution
    inv_A(:, i) = R \ y; % Backward substitution
end

% Arnoldi procedure of m
for j = 1:m
    w = inv_A * V(:, j);
    for i = 1:j
        H(i, j) = transpose(V(:, i)) * w;
        w = w - (H(i, j) * V(:, i));
    end
    % Reorthogonalization of w
    for i = 1:j
        x = transpose(V(:, i)) * w;
        w = w - (x * V(:, i));
        H(i, j) = H(i, j) + x;
    end
    H(j+1, j) = norm(w);
    if H(j+1, j) == 0
        break
    end
    V(:, j+1) = w / H(j+1, j);
end

% Calculate residual
res_1 = norm((inv_A * V(:,1:m)) - (V * H));
res_2 = norm(eye(m+1) - (transpose(V) * V));

k = 1;
[ritz, lam, rel] = plot_shift_rel_error(A, H, m, k, tau);

