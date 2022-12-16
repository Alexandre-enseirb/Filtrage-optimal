clear; clc; close all; dbstop if error;
addpath("generateurs/", "filtrage");

%% Parametres
params.alpha = 1;  % croissance de la fonction d'auto-corrélation
params.N = 500;    % nombre d'observations
params.T = 1;      % période d'echantillonage 

% variance de l'acceleration en x et en y a l'instant 0
params.sigma_2m = struct("x", 1.2, "y", 0.9);
nb_var_etat = 3;

% bruit de mesure
sigma2_measure = 200;

% Bias des paramètres du modèle au filtrage: de 20% à 180%
sigma2_errors = linspace(0.2*sigma2_measure, 1.8*sigma2_measure, 100);
alpha_errors  = linspace(0.2*params.alpha, 1.8*params.alpha, 100);

%% Une realisation
% simulation d'un mouvement singer
X.x = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.x);
X.y = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.y);

% ajout du bruit de mesure
R = get_R(sigma2_measure, nb_var_etat);  % cov bruit de mesure
Y.x = X.x + R * rand(size(X.x));
Y.y = X.y + R * rand(size(X.x));

% paramètres simulation ideaux
phi = get_phi(params.alpha, params.T);
H   = eye(nb_var_etat);
Qx  = get_Q(params.sigma_2m.x, params.alpha, params.T);  % cov bruit modele
Qy  = get_Q(params.sigma_2m.y, params.alpha, params.T);  % cov bruit modele

% filtrage ideal
[X_hat_ideal.x, P_ideal.x] = kalman(Y.x, phi, H, Qx, R);
[X_hat_ideal.y, P_ideal.y] = kalman(Y.y, phi, H, Qx, R);

%% Erreurs sur les paramètres
% la matrice de covariance du bruit de mesure
Nb_err.measure = length(sigma2_errors);
errors.measure = zeros(3, Nb_err.measure);

for i = 1:Nb_err.measure
    R_error = get_R(sigma2_errors(i), nb_var_etat);
    [X_hat.x, P.x] = kalman(Y.x, phi, H, Qx, R_error);
    errors.measure(:, i) = vecnorm(X_hat_ideal.x - X_hat.x, 2, 2);
end

% la matrice de covariance du bruit de modele
Nb_err.model = length(alpha_errors);
errors.model = zeros(3, Nb_err.model);

for i = 1:Nb_err.model
    Qx_error = get_Q(params.sigma_2m.x, alpha_errors(i), params.T);
    [X_hat.x, P.x] = kalman(Y.x, phi, H, Qx_error, R);
    errors.model(:, i) = vecnorm(X_hat_ideal.x - X_hat.x, 2, 2);
end

% les matrices de covariance du bruit de mesure et de modele
errors.cross = zeros(3, Nb_err.measure, Nb_err.model);
for i = 1:Nb_err.measure
    R_error = get_R(sigma2_errors(i), nb_var_etat);
    for j = 1:Nb_err.model
        Qx_error = get_Q(params.sigma_2m.x, alpha_errors(j), params.T);
        [X_hat.x, P.x] = kalman(Y.x, phi, H, Qx_error, R_error);
        errors.cross(:, i, j) = vecnorm(X_hat_ideal.x - X_hat.x, 2, 2);
    end
end

%% affichage
close all; clc
subName = ["Position", "Vitesse", "Acceleration"];
colors  = ["#127BCA", "#D95319", "#FF0000"];

% bruit de mesure
figure("Name", "Erreur quadratique au filrage ideal - Bruit de mesure", "Position", get(0, "ScreenSize"));
for i = 1:3
    subplot(1, 3, i)
    hold all
    plot(sigma2_errors, errors.measure(i, :), Color=colors(1))
    plot([sigma2_measure, sigma2_measure], [min(errors.measure(i, :)), max(errors.measure(i, :))], Color=colors(2))
    grid on
    xlabel("variance du bruit de mesure")
    legend(["Erreur", "frontiere"], "Location", "northeast")
    axis("square")
    title(subName(i))
end

% Bruit de modele
figure("Name", "Erreur quadratique au filrage ideal - Bruit de modele", "Position", get(0, "ScreenSize"))
for i = 1:3
    subplot(1, 3, i)
    hold all
    plot(alpha_errors, errors.model(i, :), Color=colors(1))
    plot([params.alpha, params.alpha], [min(errors.model(i, :)), max(errors.model(i, :))], Color=colors(2))
    grid on
    xlabel("alpha")
    legend(["Erreur", "frontiere"], "Location", "northeast")
    axis("square")
    title(subName(i))
end

% Bruit de mesure et de modele
figure("Name", "Erreurs croisées", "Position", get(0, "ScreenSize"))
for i = 1:3
    subplot(2, 2, i)
    hold all
    imagesc(sigma2_errors, alpha_errors, squeeze(errors.cross(i, :, :)))
    plot([sigma2_measure, sigma2_measure], [alpha_errors(1), alpha_errors(end)], Color=colors(3))
    plot([sigma2_errors(1), sigma2_errors(end)], [params.alpha, params.alpha], Color=colors(3))
    title(subName(i))
    xlabel("variance du bruit de mesure")
    ylabel("alpha")
    axis("square", "tight")
    colorbar
end
