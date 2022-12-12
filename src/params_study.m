clear; clc; close all; dbstop if error;

addpath("generateurs/", "filtrage");

%% Parametres
params.alpha = 1;  % croissance de la fonction d'auto-corrélation
params.N = 500;    % nombre d'observations
params.T = 1;      % période d'echantillonage
% variance de l'acceleration en x et en y a l'instant 0
params.sigma_2m = struct("x", 1.2, "y", 0.9);
nb_var_etat = 3;

% Bias des paramètres du modèle au filtrage


%% Une realisation
% simulation d'un mouvement singer
X.x = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.x);
X.y = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.y);

% ajout du bruit d'observation
R = 200*eye(nb_var_etat);
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
X_smooth.x = lissage(Y.x, X_hat.x, P.x, phi);
X_smooth.y = lissage(Y.y, X_hat.y, P.y, phi);