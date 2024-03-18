% Compute relative error for shirt-and-invert spectral transformation
function [ritz, lam, rel] = plot_shift_rel_error(A, H, m, k, tau)
    [n, ~] = size(A);
    rel = zeros(m, k);

    lam = eig(full(A));
    for j = 1:m
        lam(j, 1) = 1 / (lam(j, 1) - tau);
    end
    lam = sort(abs(lam), 'descend');

    for j = 1:m
        ritz = eig(H(1:j, 1:j));
        for i = 1:j
            ritz(i, 1) = 1 / ritz(i, 1);
            ritz(i, 1) = ritz(i, 1) + tau;
        end
        ritz = sort(abs(ritz), 'descend');

        for i = 1:k
            if j >= i
                rel(j, i) = rel(j, i) + abs(lam(i) - ritz(i)) / abs(lam(i));
            end
        end
    end
    
    plot(1:m, rel(1:m, 1), 'b-')
    %plot(1:m, rel(1:m, 1), 'r-', 2:m, rel(2:m, 2), 'g-', 3:m, rel(3:m, 3), 'b-');
end