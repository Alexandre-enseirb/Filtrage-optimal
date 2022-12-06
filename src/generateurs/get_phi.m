function phi = get_phi(alpha, T)
% PHI renvoie la matrice de passage pour un mouvement de Singer
    at = alpha * T;
    phi = [1  T (exp(-at) + at - 1)/alpha^2 ;...
           0  1 (1 - exp(-at))/alpha        ;...
           0  0 exp(-at)                   ];
end
