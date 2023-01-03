function [X_hat, P] = kalman(Y, phi, H, Q, R, G)
%KALMAN applique le filtre de Kalman sur le vecteur d'état observé Y,
%connaissant la matrice de transition phi, la matrice de covariance du
%bruit de modèle, le filtre H, et la matrice de covariance du bruit des observations R.
    [nb_var_etat, N] = size(Y);     % nombre d'échantillons à traiter
    X_hat = zeros(nb_var_etat, N);  % estimation des vecteurs d'état 
    X_k_apriori = zeros(nb_var_etat, 1);   % avant la connaissance de y_k
    X_k_aposter = zeros(nb_var_etat, 1);   % après la connaissance de y_k
    K_k = zeros(nb_var_etat, nb_var_etat); % gain de Kalman
    P.apriori = zeros(nb_var_etat, nb_var_etat, N);  % (3x3xN)
    P.aposter = zeros(nb_var_etat, nb_var_etat, N);  % (3x3xN)

    for k = 2:N
        X_k_apriori = phi * X_k_aposter;
        P.apriori(:, :, k) = phi * P.aposter(:, :, k-1) * phi.' + G(k)*Q*G(k).';
        K_k = P.apriori(:, :, k) * H.' / (H * P.apriori(:, :, k) * H.' + R);
        X_k_aposter = X_k_apriori + K_k * (Y(:, k) - H*X_k_apriori);
        P.aposter(:, :, k) = (eye(nb_var_etat) - K_k * H) * P.apriori(:, :, k);
        X_hat(:, k) = X_k_aposter;
    end
end

