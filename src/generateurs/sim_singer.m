function [X] = sim_singer(N, a, T, sigma2)
%SIM_SINGER simule sur N points l'evolution d'un processus a deux
%dimensions x et y selon les parametres a et T

if ~exist('sigma2', 'var')
    sigma2 = 0.0;
end

X = zeros(3,N);
mu = [0; 0; 0];

phi = [1  T (exp(-a*T) + a * T - 1)/a^2;...
    0  1      (1 - exp(-a*T))/a;     ...
    0  0          exp(-a*T)         ];

% mise a jour
for i=1:N-1
    X_tmp      = phi * X(:,i);
    Q          = chol(get_Q(sigma2, a, T));
    W          = mu + Q.'*randn(size(Q,2),1);
    X(:,i+1) = X_tmp+W;
end



