function [X] = sim_singer(N, alpha, T, sigma2, G)
%SIM_SINGER simule sur N points l'évolution d'un processus a deux
%dimensions x et y selon les parametres a et T

if ~exist('sigma2', 'var')
    sigma2 = 0.0;
end

nb_var_etat = 3;  % position, vitesse et accéleration
X = zeros(nb_var_etat, N);

phi = get_phi(alpha, T);
Q = get_Q(sigma2, alpha, T);
A = chol(Q).';  % Q = A*(A^t)

% mise a jour
for i=1:N-1
    X_next = phi * X(:, i);
    W      = A*randn(nb_var_etat, 1);
    X(:, i+1) = X_next + G(k)*W;
end
