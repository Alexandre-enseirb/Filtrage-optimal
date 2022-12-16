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

params.monte_carlo_nSim = 1;

%% MONTE CARLO

for i=1:params.monte_carlo_nSim
    errors = monte_carlo(params, sigma2_measure, params.monte_carlo_nSim);
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
