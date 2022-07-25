function[S] = special_fun2(varargin)

% Estimation of \sum_{|j| > M} (a^(j/2))/(a^{j+2} + a^2)

a = varargin{1};

if nargin > 1
    M = varargin{2};
else
    M = -1;
end

if M <= -1
    S(1) = 1/(2*a^2) + 2/(a^(5/2)+a^(3/2));
    S(2) = S(1) + 2/(a^3+a);
    j = 2;
elseif M == 0
    S(1) = 2/(a^(5/2)+a^(3/2));
    S(2) = S(1) + 2/(a^3+a);
    j = 2;
else
    S(M) = 0;
    S(M+1) = 2/(a^(2+(M+1)/2) + a^(2-(M+1)/2));
    j = M+1;
end

while abs(S(j)/S(j-1)) > 1 + 10^(-3)
    j = j+1;
    S(j) = S(j-1) + 2/(a^(2+j/2) + a^(2-j/2));
end

% Adding the error estimate (see Remark 23)
S = S(j) + 4*atan(a^(j/2))/(a^2*abs(log(a)));