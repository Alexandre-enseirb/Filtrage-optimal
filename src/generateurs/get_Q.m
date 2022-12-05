function Q = get_Q(sigma, a, T)
    at = a*T;
    q11 = 1/(2*a^5) * (2*at - 2*at^2 + 2/3*at^3 - 4*at*exp(-at) - exp(-2*at) + 1);
    q12 = 1/(2*a^4) * (at^2 + 1 + exp(-2*at) + exp(-at)*(-2 + 2*at) - 2*at);
    q13 = 1/(2*a^3) * (1 - 2*at*exp(-at) - exp(-2*at));
    q22 = 1/(2*a^3) * (2*at - 3 + 4*exp(-at) - exp(-2*at));
    q23 = 1/(2*a^2) * (1 - exp(-at))^2;
    q33 = 1/(2*a^1) * (1 - exp(-2*at));

    Q = 2*a*sigma * [q11 q12 q13;
                     q12 q22 q23;
                     q13 q23 q33];
end

