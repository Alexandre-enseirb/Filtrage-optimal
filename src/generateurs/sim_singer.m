function [X] = sim_singer(N, a, T, sigma2)
%SIM_SINGER simule sur N points l'evolution d'un processus a deux
%dimensions x et y selon les parametres a et T

if ~exist('sigma2', 'var')
    sigma2=0;
end

X = zeros(6,N);
mu = [0; 0; 0];

for i=1:N-1
    X_tmp = singer(X(:,i), a, T);
    Q     = get_Q(sigma2, a, T);
    R     = chol(Q);
    W     = mu + [randn(size(X))*R; randn(size(X))*R];
    X(:,i+1) = X_tmp+W;
end