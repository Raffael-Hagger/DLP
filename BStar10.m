function[B] = BStar10(a)

% Computation of B_n^* for n = 10 (see (4.25))

j = 2:10;
F = a.^(j/2)./(a-a.^j);
B = sum(F) + log(tanh(abs(9*log(a))/4))/(sqrt(a)*log(a));