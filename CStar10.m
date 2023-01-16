function[C] = CStar10(a)

% Computation of C_n^* for n = 10 (see (4.58))

j = 1:10;
F = 1./(a.^(j/2)+a.^(-j/2));
C = 1/(2*a^2) + 2*sum(F)/a^2 + 4*atan(a^5)/(a^2*abs(log(a)));