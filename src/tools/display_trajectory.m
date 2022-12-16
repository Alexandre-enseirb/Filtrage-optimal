function [] = display_trajectory(X,Y,X_hat,X_smooth,fig_title)
%DISPLAY_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

figure("Position", get(0, "ScreenSize"))
if exist('fig_title','var')
    sgtitle(fig_title)
end
subplot(221)
plot(X.x(1, :), X.y(1, :))
title("Trajectoire réelle")
xlabel("X")
ylabel("Y")
grid

subplot(222)
plot(Y.x(1, :), Y.y(1, :), ".")
title("Trajectoire bruitée")
xlabel("X")
ylabel("Y")
grid

subplot(223)
plot(X_hat.x(1, :), X_hat.y(1, :));
title("Trajectoire estimée")
xlabel("X")
ylabel("Y")
grid

subplot(224)
plot(X_smooth.x(1, :), X_smooth.y(1, :));
title("Trajectoire lissée")
xlabel("X")
ylabel("Y")
grid
end

