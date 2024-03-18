% Compute relative error
function rel = plot_rel_error(A, H, m, k)
    lam = eig(full(A));
    lam = sort(abs(lam), 'descend'); % Sort eigenvalues of A
    rel = zeros(m, k);

    for j = 1:m
        ritz = eig(H(1:j, 1:j));
        ritz = sort(abs(ritz), 'descend'); % Sort eigenvalues of H_j
        for i = 1:k
            if j >= i
                rel(j, i) = rel(j, i) + abs(lam(i) - ritz(i)) / abs(lam(i));
            end
        end
    end

    %plot(1:m, rel(1:m, 1), 'r-', 2:m, rel(2:m, 2), 'g-', 3:m, rel(3:m, 3), 'b-');
end