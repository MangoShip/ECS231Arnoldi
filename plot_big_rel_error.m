function rel = plot_big_rel_error(A, V, H, m, k)
    rel = zeros(m, k);
    norm_A = norm(full(A));

    for j = 1:m
        [ritz_V, ritz_D] = eig(H(1:j, 1:j));
        ritz = eig(H(1:j, 1:j));
        ritz = sort(ritz, 'descend'); % Sort eigenvalues of H_j
        for i = 1:k
            if j >= i
                [~, ind_y] = find(ritz_D==ritz(i));
                eig_vec = ritz_V(:, ind_y);
                ritz_vec = V(:,1:j) * eig_vec;
                
                rel(j, i) = rel(j, i) + (norm((A * ritz_vec) - (ritz(i) * ritz_vec)) / ((norm_A + abs(ritz(i))) * norm(ritz_vec)));
            end
        end
    end
    plot(1:m, rel(1:m, 1), 'r-', 2:m, rel(2:m, 2), 'g-', 3:m, rel(3:m, 3), 'b-');
end