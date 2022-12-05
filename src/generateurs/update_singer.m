function [X_next] = update_singer(X_prev, a, T)
%SINGER met a jour le vecteur d'etat par la methode de singer
%sans appliquer le bruit

phi = [1  T (exp(-a*T) + a * T - 1)/a^2;...
       0  1      (1 - exp(-a*T))/a;     ...
       0  0          exp(-a*T)         ];

phi_ext = [phi zeros(size(phi));...
           zeros(size(phi)) phi];

X_next = phi_ext * X_prev;

       