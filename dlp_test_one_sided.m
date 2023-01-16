function[L_c,R_c,Rho] = dlp_test_one_sided(varargin)

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
t = ((1:m) - 1/2)*pi/m;

% The estimates for \|\tilde{K}_t\|_{c,0} and \|\tilde{K}_t\|_{0,c},
% respectively (Prop. 4.8)
I_c = nimg + nimg1/lga;
F_0 = ng + ng1/lga;
F_c = ngc + ngc1/lga;
G_0 = ng1 + ng2/lga;
G_c = ngc1 + ngc2/lga;
if c > acos(a)/lga || I_c >= 1
    warning('Condition (4.21) violated.')
end
C_5 = (1 + F_c^2)^(1/4)/(pi*(1 - I_c^2))*(G_c*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + lga*(F_c+F_0)*BStar10(a));
C_6 = (1 + F_0^2)^(1/4)/(pi*(1 - I_c^2)^(5/4))*(G_c*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + 2*lga*F_c*BStar10(a));
C_4 = 2*exp(2*pi*c)*C_5*C_6;

% Estimate for \|\tilde{K}_{t,N}\|_{\infty} (see (4.39))
C_3 = (1 + F_0^2)^(1/4)/pi*(G_0*(1 + a^(1/2) + a^(-1/2)/(4*a^2)) + 2*lga*F_0*BStar10(a));

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

% Computation of C_1(M) and C_2 := \|B_N^M\|_{\infty} (see (4.34) and
% (4.35))
C_1 = (2*lga*F_0/pi)*(1 + F_0^2)^(1/4)*log(tanh(abs((M-1)*log(a))/4))/(sqrt(a)*log(a));
C_2 = norm(sum(J.*abs(Pos) + J.*abs(Neg),3),'inf');

% Summing the terms up
rho = zeros(m,1);
nu = zeros(m,2);
nu_ = zeros(m-1,2);
progress = waitbar(0,'Progress');
for k = 1:m
    A = Z + sum(exp(1i*J*t(k)).*Pos + exp(-1i*J*t(k)).*Neg,3);

% Computing the spectral radii of A_{t_k,N}^M for k = 1,...,m
    rho(k) = abs(eigs(A,1));

% Computing the the smallest singular value of (A_{t_k,N}^M - \mu_{k,l} I),
% where \mu_{k,l} is chosen as in Thm. 4.15
    phi = 0;
    nu(k,1) = 1/norm(inv(A - rho_0*eye(N)),'inf');
    l = 1;
    while phi < 2*pi
        phi = phi + nu(k,l)/(2*rho_0);
        mu = rho_0*exp(1i*phi);
        nu(k,l+1) = 1/norm(inv(A - mu*eye(N)),'inf');
        nu_(k,l) = nu(k,l)/4 + nu(k,l+1)/2;
        l = l+1;
    end
    text = strcat('Progress:',32,'N',32,'=',32,sprintf('%1.0f',N),44,32,'k',32,'=',32,sprintf('%1.0f',k));
    waitbar(k/(m+1),progress,text)
end
delete(progress);

% Check if \rho(A_{t_k,N}^M) < \rho_0 for all k = 1,...,m
Rho = max(rho);
if Rho >= rho_0
    warning('Spectral radius condition violated')
end

% Computation of L_c and R_c in (4.43)
L_c = C_4/(exp(2*pi*c*N) - 1);
R = min(nu_(nu_>0));
R_c = rho_0^2*(1 + C_3*(R - C_1 - C_2*pi/(2*m))^(-1))^(-1);

% Warning if M, m or n are chosen too small so that R_c is
% negative
if R_c < 0
    warning('R_c is negative (N = %d)',N)
end