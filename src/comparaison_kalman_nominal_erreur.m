clear; clc; close all; dbstop if error;

addpath("generateurs/", "filtrage", "tools");

%% Parametres
params.alpha.nominal = 1;  % croissance de la fonction d'auto-corrélation
params.N = 500;            % nombre d'observations
params.T = 1;              % période d'echantillonage

params.alpha.sous_estime = 0.2 * params.alpha.nominal;
params.alpha.sur_estime  = 1.8 * params.alpha.nominal;
G = 1*ones(params.N,1);

% variance de l'acceleration en x et en y a l'instant 0
params.sigma_2m = struct("x", 1.2, "y", 0.9);

nb_var_etat = 3;

% bruit de mesure
sigma2_measure.nominal     = 200;
sigma2_measure.sous_estime = 0.2 * sigma2_measure.nominal;
sigma2_measure.sur_estime  = 1.8 * sigma2_measure.nominal;

nb_combinaisons = 7; % combinaisons de parametres

%% Une realisation
% simulation d'un mouvement singer
X.x = sim_singer(params.N, params.alpha.nominal, params.T, params.sigma_2m.x, G);
X.y = sim_singer(params.N, params.alpha.nominal, params.T, params.sigma_2m.y, G);

% ajout du bruit d'observation
R.nominal = get_R(sigma2_measure.nominal, nb_var_etat);
Y.x = X.x + R.nominal * rand(size(X.x));
Y.y = X.y + R.nominal * rand(size(X.x));

R.sous_estime = get_R(sigma2_measure.sous_estime, nb_var_etat);
R.sur_estime  = get_R(sigma2_measure.sur_estime, nb_var_etat);

% paramètre simulation
phi = get_phi(params.alpha.nominal, params.T);
H   = eye(nb_var_etat);
Qx.nominal  = get_Q(params.sigma_2m.x, params.alpha.nominal, params.T);
Qy.nominal  = get_Q(params.sigma_2m.y, params.alpha.nominal, params.T);

Qx.sous_estime = get_Q(params.sigma_2m.x, params.alpha.sous_estime, params.T);
Qy.sous_estime = get_Q(params.sigma_2m.y, params.alpha.sous_estime, params.T);

Qx.sur_estime = get_Q(params.sigma_2m.x, params.alpha.sur_estime, params.T);
Qy.sur_estime = get_Q(params.sigma_2m.y, params.alpha.sur_estime, params.T);

%% -- filtrage
Qx_cell = {Qx.nominal, Qx.sous_estime, Qx.nominal, Qx.sous_estime, Qx.sur_estime, Qx.nominal, Qx.sur_estime};
Qy_cell = {Qy.nominal, Qy.sous_estime, Qy.nominal, Qy.sous_estime, Qy.sur_estime, Qy.nominal, Qy.sur_estime};
R_cell  = {R.nominal, R.nominal, R.sous_estime, R.sous_estime, R.nominal, R.sur_estime, R.sur_estime};

for i=1:nb_combinaisons
    % -- extraction des parametres
    R_i  = cell2mat(R_cell(i));
    Qx_i = cell2mat(Qx_cell(i));
    Qy_i = cell2mat(Qy_cell(i));

    % -- filtrage
    [X_hat(i).x, P(i).x] = kalman(Y.x, phi, H, Qx_i, R_i, G);
    [X_hat(i).y, P(i).y] = kalman(Y.y, phi, H, Qy_i, R_i, G);

    % -- lissage
    X_smooth(i).x = lissage(X_hat(i).x, P(i).x, phi);
    X_smooth(i).y = lissage(X_hat(i).y, P(i).y, phi);
end

%% -- Affichage
titles = ["Cas nominal"...
    "Sous-estimation de alpha"...
    "Sous-estimation de R"...
    "Sous-estimation de alpha et R"...
    "Sur-estimation de alpha"...
    "Sur-estimation de R"...
    "Sur-estimation de alpha et R"];

for i=1:nb_combinaisons
        display_trajectory(X,Y,X_hat(i),X_smooth(i),titles(i));
end