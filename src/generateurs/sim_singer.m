function [X] = sim_singer(N, a, T, sigma2)
    %SIM_SINGER simule sur N points l'evolution d'un processus a deux
    %dimensions x et y selon les parametres a et T
    
    if ~exist('sigma2', 'var')
        sigma2 = struct("x", 0.0, "y", 0.0);
    end
    
    X = zeros(6,N);
    mu = [0; 0; 0];
    
    for i=1:N-1
        X_tmp = update_singer(X(:,i), a, T);
        Q_x   = get_Q(sigma2.x, a, T);
        Q_y   = get_Q(sigma2.y, a, T);
        R     = [chol(Q_x);chol(Q_y)];
        W     = repmat(mu,2,1) + [R*randn(size(R,2),1)];% R*randn(size(R,2),1)];
        X(:,i+1) = X_tmp+W;
    end