function X_hat = kalman(Y, phi, H, Q, R)
%KALMAN applique le filtre de Kalman sur le vecteur d'état observé Y,
%connaissant la matrice de transition phi, la matrice de covariance du
%bruit de modèle, le filtre H, et la matrice de covariance du bruit des observations R.
    [nb_var_etat, N] = size(Y);     % nombre d'échantillons à traiter
    X_hat = zeros(nb_var_etat, N);  % estimation des vecteurs d'état 
    X_k_apriori = zeros(nb_var_etat, 1);   % avant la connaissance de y_k
    X_k_aposter = zeros(nb_var_etat, 1);   % après la connaissance de y_k
    K_k = zeros(nb_var_etat, nb_var_etat); % gain de Kalman
    P_apriori = zeros(nb_var_etat, nb_var_etat);  % (3x3)
    P_aposter = zeros(nb_var_etat, nb_var_etat);  % (3x3)

    for k = 2:N
        X_k_apriori = phi * X_k_aposter;
        P_apriori = phi * P_aposter * phi.' + Q*Q.';  % Q*Q' <-> Q*G*Q' ???
        K_k = P_apriori * H.' / (H * P_apriori * H.' + R);
        X_k_aposter = X_k_apriori + K_k * (Y(:, k) - H*X_k_apriori);
        P_aposter = (eye(nb_var_etat) - K_k * H) * P_apriori;
        X_hat(:, k) = X_k_aposter;
    end
end

