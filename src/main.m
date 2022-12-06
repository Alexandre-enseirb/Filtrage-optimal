clear; clc; close all; dbstop if error;

addpath("generateurs/");

%% 
alpha = 0.1; % croissance de la fonction d'auto-corrélation
N = 500;     % nombre d'observations
T = 1;       % période d'échantillonage
sigma_2m = struct("x", 1.2, "y", 0.9);


%%
mu = [0, 0, 0];
%Q = get_Q(sigma_2m.x, alpha, T);
%R = chol(Q);
%z = repmat(mu, N, 1) + randn(N, 3)*R;

X = sim_singer(N, alpha, T, sigma_2m);

X_pos_x = X(1,:);
X_pos_y = X(4,:);

figure
plot(X_pos_x, X_pos_y, 'Color', '#4DBEEE', 'LineStyle','--');
xlabel("X position");
ylabel("Y position");
legend("Singer movement");
title("Simulation of a Singer movement over 500 iterations.");