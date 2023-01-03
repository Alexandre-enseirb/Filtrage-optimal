clear; clc; close all; dbstop if error;

addpath("generateurs/", "filtrage", "tools");

%% PARAMETRES

params.alpha = 1;  % croissance de la fonction d'auto-corrélation
params.N = 500;    % nombre d'observations
params.T = 1;      % période d'echantillonage

% variance de l'acceleration en x et en y a l'instant 0
params.sigma_2m = struct("x", 1.2, "y", 0.9);
nb_var_etat = 3;
G = 1*ones(1,params.N);
% bruit de mesure
sigma2_measure = 200;

nSim = 300;

%% MONTE-CARLO

dist_bruite = zeros(1,params.N);
dist_kalman = zeros(1,params.N);
dist_lissee = zeros(1,params.N);

for i=1:nSim
    % simulation d'un mouvement singer
    X.x = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.x, G);
    X.y = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.y, G);
    
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
    [X_hat.x, P.x] = kalman(Y.x, phi, H, Qx, R, G);
    [X_hat.y, P.y] = kalman(Y.y, phi, H, Qx, R, G);
    
    % lissage
    X_smooth.x = lissage(X_hat.x, P.x, phi);
    X_smooth.y = lissage(X_hat.y, P.y, phi);

    % erreur quadratique
    dist_bruite = dist_bruite + sqrt( (X.x(1,:) - Y.x(1,:)).^2        + (X.y(1,:) - Y.y(1,:)).^2 );
    dist_kalman = dist_kalman + sqrt( (X.x(1,:) - X_hat.x(1,:)).^2    + (X.y(1,:) - X_hat.y(1,:)).^2 );
    dist_lissee = dist_lissee + sqrt( (X.x(1,:) - X_smooth.x(1,:)).^2 + (X.y(1,:) - X_smooth.y(1,:)).^2 );
end

dist_bruite_moy = sum(dist_bruite)/(nSim*params.N);
dist_kalman_moy = sum(dist_kalman)/(nSim*params.N);
dist_lissee_moy = sum(dist_lissee)/(nSim*params.N);

dist_bruite_pts = dist_bruite/nSim;
dist_kalman_pts = dist_kalman/nSim;
dist_lissee_pts = dist_lissee/nSim;

figure(1);
% -- affichage des trajectoires
subplot(121)
plot(X.x(1,:), X.y(1,:));
xlabel("X");
ylabel("Y");
hold on;
grid on;
plot(Y.x(1,:), Y.y(1,:), '+');
title("Trajectoire réalisée et mesure bruitée");
legend(["Trajectoire de Singer", "Mesure bruitée"]);
% -- distance quadratique associee
subplot(122)
plot(1:params.N, dist_bruite_pts);
hold on;
grid on;
xlabel("K")
ylabel("{||.||^2}");
plot(1:params.N, dist_bruite_moy*ones(1,params.N));
legend(["Distance moyenne par point", "Distance moyenne totale"]);
title("Distance moyenne à la trajectoire lors de la mesure");

figure(2);
% -- affichage des trajectoires
subplot(121)
plot(X.x(1,:), X.y(1,:));
xlabel("X");
ylabel("Y");
hold on;
grid on;
plot(X_hat.x(1,:), X_hat.y(1,:));
title("Trajectoire réalisée et filtrage de Kalman");
legend(["Trajectoire de Singer", "Trajectoire estimée"]);
% -- distance quadratique associee
subplot(122)
plot(1:params.N, dist_kalman_pts);
hold on;
grid on;
xlabel("K")
ylabel("{||.||^2}");
plot(1:params.N, dist_kalman_moy*ones(1,params.N));
plot(1:params.N, 1.02*dist_kalman_moy*ones(1,params.N), 'r--');
plot(1:params.N, 0.98*dist_kalman_moy*ones(1,params.N), 'r--');
legend(["Distance moyenne par point", "Distance moyenne totale", "bornes de l'intervalle de confiance"]);
title("Distance moyenne à la trajectoire lors de l'estimation");

figure(3);
% -- affichage des trajectoires
subplot(121)
plot(X.x(1,:), X.y(1,:));
xlabel("X");
ylabel("Y");
hold on;
grid on;
plot(X_smooth.x(1,:), X_smooth.y(1,:));
title("Trajectoire réalisée et lissage de Kalman");
legend(["Trajectoire de Singer", "Trajectoire lissée"]);
% -- distance quadratique associee
subplot(122)
plot(1:params.N, dist_lissee_pts);
hold on;
grid on;
xlabel("K")
ylabel("{||.||^2}");
plot(1:params.N, dist_lissee_moy*ones(1,params.N));
legend(["Distance moyenne par point", "Distance moyenne totale"]);
title("Distance moyenne à la trajectoire lors du lissage");