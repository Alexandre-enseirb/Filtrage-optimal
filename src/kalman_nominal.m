clear; clc; close all; dbstop if error;

addpath("generateurs/", "filtrage", "tools");

%% Parametres
params.alpha = 1;  % croissance de la fonction d'auto-corrélation
params.N = 500;    % nombre d'observations
params.T = 1;      % période d'echantillonage

% variance de l'acceleration en x et en y a l'instant 0
params.sigma_2m = struct("x", 1.2, "y", 0.9);
nb_var_etat = 3;

% bruit de mesure
sigma2_measure = 200;

%% Une realisation
% simulation d'un mouvement singer
X.x = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.x);
X.y = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.y);

% ajout du bruit d'observation
R = get_R(sigma2_measure, nb_var_etat);
Y.x = X.x + R * rand(size(X.x));
Y.y = X.y + R * rand(size(X.x));

% paramètre simulation
phi = get_phi(params.alpha, params.T);
H   = eye(nb_var_etat);
Qx  = get_Q(params.sigma_2m.x, params.alpha, params.T);
Qy  = get_Q(params.sigma_2m.y, params.alpha, params.T);

% filtrage
[X_hat.x, P.x] = kalman(Y.x, phi, H, Qx, R);
[X_hat.y, P.y] = kalman(Y.y, phi, H, Qx, R);

% lissage
X_smooth.x = lissage(X_hat.x, P.x, phi);
X_smooth.y = lissage(X_hat.y, P.y, phi);

%% -- Affichage
display_trajectory(X,Y,X_hat,X_smooth,"Filtrage de Kalman nominal");

