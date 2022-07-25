function[up_est,low_est,Rho] = dlp_test_one_sided(varargin)

% Setting some default options
if nargin > 0
    N = varargin{1};
else
    N = 256;
end
if nargin > 1
    a = varargin{2}; % a := \alpha
else
    a = 1/2;
end
if nargin > 2
    c = varargin{3};
else
    c = 0.03;
end
if nargin > 3
    m = varargin{4};
else
    m = 400;
end
if nargin > 4
    M = varargin{5};
else
    M = 25;
end
if nargin > 5
    rho_0 = varargin{6};
else
    rho_0 = 1/2;
end
if nargin > 6
    g = varargin{7};
    ng = varargin{8}; % sup norm of g
    ng1 = varargin{9}; % sup norm of g'
    ng2 = varargin{10}; % sup norm of g''
    nimg = varargin{11}; % c norm of Im(g)
    nimg1 = varargin{12}; % c norm of Im(g')
    ngc = varargin{13}; % c norm of g
    ngc1 = varargin{14}; % c norm of g'
    ngc2 = varargin{15}; % c norm of g''
else
    g = @(x) sin(pi*x).^2;
    ng = 1;
    ng1 = pi;
    ng2 = 2*pi^2;
    nimg = (1/2)*sinh(2*pi*c);
    nimg1 = pi*sinh(2*pi*c);
    ngc = cosh(pi*c)^2;
    ngc1 = pi*cosh(2*pi*c);
    ngc2 = 2*pi^2*cosh(2*pi*c);
end

if c > acos(a)/abs(log(a)) || nimg + nimg1/abs(log(a)) >= 1
    warning('Conditions on c violated.')
end

% Derivatives of g
syms x;
g1 = matlabFunction(diff(g(x),x),'Vars',x);
g2 = matlabFunction(diff(g1(x),x),'Vars',x);
clear x;

% The function generating the graph and its derivatives
f = @(x) x.*g(log(x)/log(a));
f1 = @(x) g(log(x)/log(a)) + g1(log(x)/log(a))/log(a);
f2 = @(x) g1(log(x)/log(a))./(x*log(a)) + g2(log(x)/log(a))./(x*(log(a))^2);

% Abrreviate abs(log(a)) by lga for readability
lga = abs(log(a));

% Kernel of the DLP operator after being transformed to an operator on \R.
% kernel_sing is the kernel on the diagonal.
kernel = @(x,y) (1/(2*pi))*(((a.^y-a.^x).*f1(a.^y) + f(a.^x) - f(a.^y))./((a.^x-a.^y).^2 + (f(a.^x)-f(a.^y)).^2)).*((1+f1(a.^x).^2)./(1+f1(a.^y).^2)).^(1/4).*a.^((x+y)/2)*lga;
kernel_sing = @(x) (1/(2*pi))*f2(a.^x)./(2+2*f1(a.^x)^2).*a.^x*lga;

% Fix t_k
t = (0:m)*pi/m;

% The estimates for \|\tilde{K}_t\|_{c,0} and \|\tilde{K}_t\|_{0,c},
% respectively
K_1 = (1 + (ngc + ngc1/lga)^2)^(1/4)/(pi*(1 - (nimg + nimg1/lga)^2))*((ngc1 + ngc2/lga)*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + (lga*(ngc + ng) + ngc1 + ng1)*special_fun1(a));
K_2 = (1 + (ng + ng1/lga)^2)^(1/4)/(pi*(1 - (nimg + nimg1/lga)^2)^(5/4))*((ngc1 + ngc2/lga)*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + 2*(lga*ngc + ngc1)*special_fun1(a));

% Estimate for \|\tilde{K}_t\|_{\infty}
K_3 = (1 + (ng + ng1/lga)^2)^(1/4)/pi*((ng1 + ng2/lga)*(1 + a^(1/2) + a^(-1/2)/(4*a^2)) + 2*(lga*ng + ng1)*special_fun1(a));

% Computation of B_{t_k,N}^M for k = 0,...,m. As the summands are
% identical for all k except for the exp(ijt_k), we first compute the
% summands and then add them up later.

% The 0-th term
x = 1/(2*N):1/N:1-1/(2*N);
y = x;
[Y,X] = meshgrid(y,x);
Z = kernel(X,Y)/N;
for p = 1:N
    Z(p,p) = kernel_sing(x(p))/N;
end

% The j-th term, j = -M,...,-1,1,...,M
j = 1:M;
[Y,X,J] = meshgrid(y,x,j);
Pos = kernel(X+J,Y)/N;
Neg = kernel(X-J,Y)/N;

% Estimates for the constants C_1(M) and C_2(M,N)
C_1 = (2/pi)*(1 + (ng + ng1/lga)^2)^(1/4)*(lga*ng + ng1)*special_fun1(a,M);
C_2 = norm(sum(J.*abs(Pos) + J.*abs(Neg),3),'inf');

% Summing the terms up
rho = zeros(m+1,1);
nu = zeros(m+1,2);
nu_ = zeros(m,2);
progress = waitbar(0,'Progress');
for k = 1:m+1
    B = Z + sum(exp(1i*J*t(k)).*Pos + exp(-1i*J*t(k)).*Neg,3);

% Computing the spectral radii of B_{t_k,N}^M for k = 0,...,m
    rho(k) = abs(eigs(B,1));

% Computing the the smallest singular value of (B_{t_k,N}^M - \lambda_l I),
% where \lambda_l is chosen adaptively
    phi = 0;
    nu(k,1) = 1/norm(inv(B - rho_0*eye(N)),'inf');
    l = 1;
    while phi < 2*pi
        phi = phi + nu(k,l)/(2*rho_0);
        lambda = rho_0*exp(1i*phi);
        nu(k,l+1) = 1/norm(inv(B - lambda*eye(N)),'inf');
        nu_(k,l) = nu(k,l)/4 + nu(k,l+1)/2;
        l = l+1;
    end
    text = strcat('Progress:',32,'N',32,'=',32,sprintf('%1.0f',N),44,32,'k',32,'=',32,sprintf('%1.0f',k));
    waitbar(k/(m+1),progress,text)
end
delete(progress);

% Check if \rho(B_{t_k,N}^M) < \rho_0 for all k = 0,...,m
Rho = max(rho);
if Rho >= rho_0
    warning('Spectral radius condition violated')
end

% Upper estimate for \|(\tilde{K}_t - \tilde{K}_{t,N})\tilde{K}_t\|_{\infty}
% (Proposition 14)
up_est = 2*exp(2*pi*c)*K_1*K_2/(exp(2*pi*c*N) - 1);

% Lower estimate for the RHS in Theorem 18
Nu = min(nu_(nu_>0));
low_est = rho_0^2*(1 + K_3*(Nu - C_2*pi/(2*m) - C_1)^(-1))^(-1);

% Warning if M,m or n are chosen too small and thus the lower estimate gets
% negative
if Nu < C_2*pi/(2*m)
    warning('m is chosen too small')
end
if Nu < C_1
    warning('M is chosen too small')
end