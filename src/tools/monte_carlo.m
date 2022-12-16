function [errors] = monte_carlo(params, sigma2_measure, nSim)

nb_var_etat = 3;
% Bias des paramètres du modèle au filtrage: de 20% à 180%
sigma2_errors = linspace(0.2*sigma2_measure, 1.8*sigma2_measure, 100);
alpha_errors  = linspace(0.2*params.alpha, 1.8*params.alpha, 100);

% paramètres simulation ideaux
phi = get_phi(params.alpha, params.T);
H   = eye(nb_var_etat);
Qx  = get_Q(params.sigma_2m.x, params.alpha, params.T);  % cov bruit modele
R = get_R(sigma2_measure, nb_var_etat);  % cov bruit de mesure

Nb_err.measure = length(sigma2_errors);
errors.measure = zeros(3, Nb_err.measure);
Nb_err.model = length(alpha_errors);
errors.model = zeros(3, Nb_err.model);
errors.cross = zeros(3, Nb_err.measure, Nb_err.model);
for s=1:nSim
    %% Une realisation
    % simulation d'un mouvement singer
    X.x = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.x);

    % ajout du bruit de mesure
    Y.x = X.x + R * rand(size(X.x));

    % filtrage ideal
    [X_hat_ideal.x, P_ideal.x] = kalman(Y.x, phi, H, Qx, R);

    %% Erreurs sur les paramètres
    % la matrice de covariance du bruit de mesure


    for i = 1:Nb_err.measure
        R_error = get_R(sigma2_errors(i), nb_var_etat);
        [X_hat.x, P.x] = kalman(Y.x, phi, H, Qx, R_error);
        errors.measure(:, i) = errors.measure(:, i) + vecnorm(X_hat_ideal.x - X_hat.x, 2, 2);
    end

    % la matrice de covariance du bruit de modele


    for i = 1:Nb_err.model
        Qx_error = get_Q(params.sigma_2m.x, alpha_errors(i), params.T);
        [X_hat.x, P.x] = kalman(Y.x, phi, H, Qx_error, R);
        errors.model(:, i) = errors.model(:, i) + vecnorm(X_hat_ideal.x - X_hat.x, 2, 2);
    end

    % les matrices de covariance du bruit de mesure et de modele

    for i = 1:Nb_err.measure
        R_error = get_R(sigma2_errors(i), nb_var_etat);
        for j = 1:Nb_err.model
            Qx_error = get_Q(params.sigma_2m.x, alpha_errors(j), params.T);
            [X_hat.x, P.x] = kalman(Y.x, phi, H, Qx_error, R_error);
            errors.cross(:, i, j) = errors.cross(:, i, j) + vecnorm(X_hat_ideal.x - X_hat.x, 2, 2);
        end
    end

end
errors.measure = errors.measure ./ nSim;
errors.model   = errors.model   ./ nSim;
errors.cross   = errors.cross   ./ nSim;
end