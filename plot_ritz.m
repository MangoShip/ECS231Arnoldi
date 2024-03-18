% Computing ritz and plotting
function plot_ritz(A, H, m)
    lam = eig(full(A));
    ritz = eig(H(1:m, :));

    plot(real(lam), imag(lam), 'r+', real(ritz), imag(ritz), 'bo');
end