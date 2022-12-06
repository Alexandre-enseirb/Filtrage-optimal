clear; clc; close all; dbstop if error;

addpath("generateurs/");

%% Parametres
params.alpha = [0.1 1 10]; % croissance de la fonction d'auto-corrélation
params.N = 500;            % nombre d'observations
params.T = 1;              % période d'échantillonage
% variance de l'acceleration en x et en y a l'instant 0
params.sigma_2m = struct("x", 1.2, "y", 0.9);


%% Une realisation

X.x = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.x);
X.y = sim_singer(params.N, params.alpha(1), params.T, params.sigma_2m.y);

X_pos_x = X.x(1,:);  % position
X_vit_x = X.x(2,:);  % vitesse
X_acc_x = X.x(3,:);  % acceleration
X_pos_y = X.y(1,:);
X_vit_y = X.y(2,:);
X_acc_y = X.y(3,:);
time_axis = linspace(1, params.N, params.N);

figure(1)
sgtitle("Représentation d'une réalisation d'un mouvement de Singer")
subplot(221)
plot(X_pos_x, X_pos_y, 'Color', '#4DBEEE', 'LineStyle','--');
xlabel("X");
ylabel("Y");
grid on;
legend("Mouvement de Singer");
title("Evolution de la position");

subplot(222)
% courbes
plot(time_axis, X_vit_x, 'Color', '#4DBEEE');
hold on;
grid on;
plot(time_axis, X_vit_y, 'Color', '#D95319');
xlabel("Axe des temps");
ylabel("vitesse");
legend(["{v_x}" "{v_y}"]);
title("Evolution de la vitesse");

subplot(212)
plot(time_axis, X_acc_x, 'Color', '#4DBEEE');
hold on;
grid on;
plot(time_axis, X_acc_y, 'Color', '#D95319');
xlabel("Axe des temps");
ylabel("accélération");
legend(["{a_x}" "{a_y}"]);
title("Evolution de l'accélération");



%% multiples realisations
nReal = 5;

X.x = zeros(3, params.N);
X.y = zeros(3, params.N);

for i=1:length(params.alpha)
    plot_realisations(nReal, params, i);
end

function [] = plot_realisations(nReal, params, idx)

figure
colors = ["#4DBEEE",...
          "#D95319",...
          "#EDB120",...
          "#0072BD",...
          "#77AC30"];

for i=1:nReal
    X.x(:,:) = sim_singer(params.N, params.alpha(idx), params.T, params.sigma_2m.x);
    X.y(:,:) = sim_singer(params.N, params.alpha(idx), params.T, params.sigma_2m.y);
    pos_x = X.x(1,:);
    pos_y = X.y(1,:);
    plot(pos_x, pos_y, 'Color', colors(i), 'LineStyle', '--');
    hold on;
end

title("Plusieurs réalisations d'un mouvement de Singer ({\alpha = }" + num2str(params.alpha(idx)) + ")");
xlabel("X");
ylabel("Y");
legend("Realisation " + num2cell(1:nReal));
grid on;

end