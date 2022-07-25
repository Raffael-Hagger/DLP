function[S] = special_fun1(varargin)

% Estimation of \sum_{j = M+1}^\infty (a^(j/2))/(a - a^j)

a = varargin{1};

if nargin > 1
    M = varargin{2};
else
    M = 1;
end

S(M) = 0;
S(M+1) = a^((M+1)/2)/(a-a^(M+1));
j = M+1;
while abs(S(j)/S(j-1)) > 1 + 10^(-3)
    j = j+1;
    S(j) = S(j-1) + a^(j/2)/(a-a^j);
end

% Adding the error estimate (see Remark 17)
S = S(j) - 2*log(tanh(abs((j-1)*log(a))/4))/(sqrt(a)*abs(log(a)));