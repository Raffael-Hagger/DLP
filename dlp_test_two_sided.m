function[up_est,low_est,Rho] = dlp_test_two_sided(varargin)

% Setting some default options
if nargin > 0
    N = varargin{1};
else
    N = 32;
end
if nargin > 1
    a = varargin{2}; % a := \alpha
else
    a = 3/4;
end
if nargin > 2
    c = varargin{3};
else
    c = 0.05;
end
if nargin > 3
    m = varargin{4};
else
    m = 400;
end
if nargin > 4
    M = varargin{5};
else
    M = 50;
end
if nargin > 5
    rho_0 = varargin{6};
else
    rho_0 = 1/2;
end
if nargin > 6
    g_p = varargin{7}; % g_+
    g_m = varargin{8}; % g_-
    ng_p = varargin{9}; % sup norm of g_+
    ng_m = varargin{10}; % sup norm of g_-
    ng1_p = varargin{11}; % sup norm of g_+'
    ng1_m = varargin{12}; % sup norm of g_-'
    ng2_p = varargin{13}; % sup norm of g_+''
    ng2_m = varargin{14}; % sup norm of g_-''
    nimg_p = varargin{15}; % c norm of Im(g_+)
    nimg_m = varargin{16}; % c norm of Im(g_-)
    nimg1_p = varargin{17}; % c norm of Im(g_+')
    nimg1_m = varargin{18}; % c norm of Im(g_-')
    ngc_p = varargin{19}; % c norm of g_+
    ngc_m = varargin{20}; % c norm of g_-
    ngc1_p = varargin{21}; % c norm of g_+'
    ngc1_m = varargin{22}; % c norm of g_-'
    ngc2_p = varargin{23}; % c norm of g_+''
    ngc2_m = varargin{24}; % c norm of g_-''
else
    g_p = @(x) 1;
    g_m = @(x) 1;
    ng_p = 1;
    ng_m = 1;
    ng1_p = 0;
    ng1_m = 0;
    ng2_p = 0;
    ng2_m = 0;
    nimg_p = 0;
    nimg_m = 0;
    nimg1_p = 0;
    nimg1_m = 0;
    ngc_p = 1;
    ngc_m = 1;
    ngc1_p = 0;
    ngc1_m = 0;
    ngc2_p = 0;
    ngc2_m = 0;
end

if c > acos(a)/abs(log(a)) || nimg_p + nimg1_p/abs(log(a)) >= 1 || nimg_m + nimg1_m/abs(log(a)) >= 1 || c*abs(log(a))*(ng_p+ngc_m) + nimg_m >= a^4 || c*abs(log(a))*(ng_m+ngc_p) + nimg_p >= a^4
    warning('Conditions on c violated.')
end

% Derivatives of g_+ and g_-
syms x;
g1_p = matlabFunction(diff(g_p(x),x),'Vars',x);
g1_m = matlabFunction(diff(g_m(x),x),'Vars',x);
g2_p = matlabFunction(diff(g1_p(x),x),'Vars',x);
g2_m = matlabFunction(diff(g1_m(x),x),'Vars',x);
clear x;

% The function generating the graph and its derivatives
f = @(x) x.*g_p(log(x)/log(a)).*(x > 0) - x.*g_m(log(-x)/log(a)).*(x < 0);
f1 = @(x) (g_p(log(x)/log(a)) + g1_p(log(x)/log(a))/log(a)).*(x > 0) - (g_m(log(-x)/log(a)) + g1_m(log(-x)/log(a))/log(a)).*(x < 0);
f2 = @(x) (g1_p(log(x)/log(a))./(x*log(a)) + g2_p(log(x)/log(a))./(x*(log(a))^2)).*(x > 0) - (g1_m(log(-x)/log(a))./(x*log(a)) + g2_m(log(-x)/log(a))./(x*(log(a))^2)).*(x < 0);

% Abrreviate abs(log(a)) by lga for readability
lga = abs(log(a));

% Kernels of the 2 x 2 DLP operator after being transformed to an operator on \R.
% kernel_sing is the kernel on the diagonal.
kernel_K_p = @(x,y) (1/(2*pi))*(((a.^y-a.^x).*f1(a.^y) + f(a.^x) - f(a.^y))./((a.^x-a.^y).^2 + (f(a.^x)-f(a.^y)).^2)).*((1+f1(a.^x).^2)./(1+f1(a.^y).^2)).^(1/4).*a.^((x+y)/2)*lga;
kernel_K_m = @(x,y) (1/(2*pi))*((-(a.^y-a.^x).*f1(-a.^y) + f(-a.^x) - f(-a.^y))./((a.^x-a.^y).^2 + (f(-a.^x)-f(-a.^y)).^2)).*((1+f1(-a.^x).^2)./(1+f1(-a.^y).^2)).^(1/4).*a.^((x+y)/2)*lga;
kernel_L_p = @(x,y) (1/(2*pi))*((-(a.^y+a.^x).*f1(-a.^y) + f(a.^x) - f(-a.^y))./((a.^x+a.^y).^2 + (f(a.^x)-f(-a.^y)).^2)).*((1+f1(a.^x).^2)./(1+f1(-a.^y).^2)).^(1/4).*a.^((x+y)/2)*lga;
kernel_L_m = @(x,y) (1/(2*pi))*(((a.^y+a.^x).*f1(a.^y) + f(-a.^x) - f(a.^y))./((a.^x+a.^y).^2 + (f(-a.^x)-f(a.^y)).^2)).*((1+f1(-a.^x).^2)./(1+f1(a.^y).^2)).^(1/4).*a.^((x+y)/2)*lga;
kernel_sing_p = @(x) (1/(2*pi))*f2(a.^x)./(2+2*f1(a.^x).^2).*a.^x*lga;
kernel_sing_m = @(x) (1/(2*pi))*f2(-a.^x)./(2+2*f1(-a.^x).^2).*a.^x*lga;

% Fix t_k
t = (0:m)*pi/m;

% The estimates for \|\tilde{K}^{\pm}_t\|_{c,0},\|\tilde{L}^{\pm}_t\|_{c,0}
% and \|\tilde{K}^{\pm}_t\|_{0,c}, \|\tilde{L}^{\pm}_t\|_{0,c},
% respectively
K_1_p = (1 + (ngc_p + ngc1_p/lga)^2)^(1/4)/(pi*(1 - (nimg_p + nimg1_p/lga)^2))*((ngc1_p + ngc2_p/lga)*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + (lga*(ngc_p + ng_p) + ngc1_p + ng1_p)*special_fun1(a));
K_1_m = (1 + (ngc_m + ngc1_m/lga)^2)^(1/4)/(pi*(1 - (nimg_m + nimg1_m/lga)^2))*((ngc1_m + ngc2_m/lga)*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + (lga*(ngc_m + ng_m) + ngc1_m + ng1_m)*special_fun1(a));
L_1_p = (2*pi)^(-1)*(ng_m*lga + ng1_m + a^(-2)*lga*max(ngc_p,ng_m))/(1 - a^(-8)*(c*lga*(ngc_p + ng_m) + nimg_p)^2)*(1 + (ngc_p + ngc1_p/lga)^2)^(1/4)*special_fun2(a);
L_1_m = (2*pi)^(-1)*(ng_p*lga + ng1_p + a^(-2)*lga*max(ngc_m,ng_p))/(1 - a^(-8)*(c*lga*(ngc_m + ng_p) + nimg_m)^2)*(1 + (ngc_m + ngc1_m/lga)^2)^(1/4)*special_fun2(a);
K_1 = [K_1_m L_1_m;L_1_p K_1_p];
K_2_p = (1 + (ng_p + ng1_p/lga)^2)^(1/4)/(pi*(1 - (nimg_p + nimg1_p/lga)^2)^(5/4))*((ngc1_p + ngc2_p/lga)*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + 2*(lga*ngc_p + ngc1_p)*special_fun1(a));
K_2_m = (1 + (ng_m + ng1_m/lga)^2)^(1/4)/(pi*(1 - (nimg_m + nimg1_m/lga)^2)^(5/4))*((ngc1_m + ngc2_m/lga)*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + 2*(lga*ngc_m + ngc1_m)*special_fun1(a));
L_2_p = (2*pi)^(-1)*(ngc_m*lga + ngc1_m + a^(-2)*lga*max(ng_p,ngc_m))/(1 - a^(-8)*(c*lga*(ng_p + ngc_m) + nimg_m)^2)*(1 + (ng_p + ng1_p/lga)^2)^(1/4)/(1 - (nimg_m + nimg1_m/lga)^2)^(1/4)*special_fun2(a);
L_2_m = (2*pi)^(-1)*(ngc_p*lga + ngc1_p + a^(-2)*lga*max(ng_m,ngc_p))/(1 - a^(-8)*(c*lga*(ng_m + ngc_p) + nimg_p)^2)*(1 + (ng_m + ng1_m/lga)^2)^(1/4)/(1 - (nimg_p + nimg1_p/lga)^2)^(1/4)*special_fun2(a);
K_2 = [K_2_m L_2_m;L_2_p K_2_p];

% Estimates for \|\tilde{K}^{\pm}_t\|_{\infty} and \|\tilde{L}^{\pm}_t\|_{\infty}
K_3_p = (1 + (ng_p + ng1_p/lga)^2)^(1/4)/pi*((ng1_p + ng2_p/lga)*(1 + a^(1/2) + a^(-1/2)/(4*a^2)) + 2*(lga*ng_p + ng1_p)*special_fun1(a));
K_3_m = (1 + (ng_m + ng1_m/lga)^2)^(1/4)/pi*((ng1_m + ng2_m/lga)*(1 + a^(1/2) + a^(-1/2)/(4*a^2)) + 2*(lga*ng_m + ng1_m)*special_fun1(a));
L_3_p = (2*pi)^(-1)*(ng_m*lga + ng1_m + a^(-2)*lga*max(ng_p,ng_m))*(1 + (ng_p + ng1_p/lga)^2)^(1/4)*special_fun2(a);
L_3_m = (2*pi)^(-1)*(ng_p*lga + ng1_p + a^(-2)*lga*max(ng_m,ng_p))*(1 + (ng_m + ng1_m/lga)^2)^(1/4)*special_fun2(a);
K_3 = [K_3_m L_3_m;L_3_p K_3_p];

% Computation of B_{t_k,N}^M for k = 0,...,m. As the summands are
% identical for all k except for the exp(ijt_k), we first compute the
% summands and then add them up later.

% The 0-th term
x = 1/(2*N):1/N:1-1/(2*N);
y = x;
[Y,X] = meshgrid(y,x);
Z_K_p = kernel_K_p(X,Y)/N;
Z_K_m = kernel_K_m(X,Y)/N;
Z_L_p = kernel_L_p(X,Y)/N;
Z_L_m = kernel_L_m(X,Y)/N;
for p = 1:N
    Z_K_p(p,p) = kernel_sing_p(x(p))/N;
    Z_K_m(p,p) = kernel_sing_m(x(p))/N;
end
Z = [Z_K_m Z_L_m;Z_L_p Z_K_p];

% The j-th term, j = -M,...,-1,1,...,M
j = 1:M;
[Y,X,J] = meshgrid(y,x,j);
Pos_K_p = kernel_K_p(X+J,Y)/N;
Pos_K_m = kernel_K_m(X+J,Y)/N;
Pos_L_p = kernel_L_p(X+J,Y)/N;
Pos_L_m = kernel_L_m(X+J,Y)/N;
Neg_K_p = kernel_K_p(X-J,Y)/N;
Neg_K_m = kernel_K_m(X-J,Y)/N;
Neg_L_p = kernel_L_p(X-J,Y)/N;
Neg_L_m = kernel_L_m(X-J,Y)/N;
Pos = [Pos_K_m Pos_L_m;Pos_L_p Pos_K_p];
Neg = [Neg_K_m Neg_L_m;Neg_L_p Neg_K_p];

% Estimates for the constants C_1(M) and C_2(M,N)
C_1_p = (1 + (ng_p + ng1_p/lga)^2)^(1/4)*(2*(lga*ng_p + ng1_p)*special_fun1(a,M) + (ng_m*lga + ng1_m + a^(-2)*lga*max(ng_p,ng_m))*special_fun2(a,M)/2);
C_1_m = (1 + (ng_m + ng1_m/lga)^2)^(1/4)*(2*(lga*ng_m + ng1_m)*special_fun1(a,M) + (ng_p*lga + ng1_p + a^(-2)*lga*max(ng_m,ng_p))*special_fun2(a,M)/2);
C_1 = max(C_1_p,C_1_m)/pi;
C_2 = norm(sum([J J;J J].*abs(Pos) + [J J;J J].*abs(Neg),3),'inf');

% Summing the terms up
rho = zeros(m+1,1);
nu = zeros(m+1,2);
nu_ = zeros(m,2);
progress = waitbar(0,'Progress');
for k = 1:m+1
    B = Z + sum(exp(1i*[J J;J J]*t(k)).*Pos + exp(-1i*[J J;J J]*t(k)).*Neg,3);

% Computing the spectral radii of B_{t_k,N}^M for k = 0,...,m
    rho(k) = abs(eigs(B,1));

% Computing the the smallest singular value of (B_{t_k,N}^M - \lambda_l I),
% where \lambda_l is chosen adaptively
    phi = 0;
    nu(k,1) = 1/norm(inv(B - rho_0*eye(2*N)),'inf');
    l = 1;
    while phi < 2*pi
        phi = phi + nu(k,l)/(2*rho_0);
        lambda = rho_0*exp(1i*phi);
        nu(k,l+1) = 1/norm(inv(B - lambda*eye(2*N)),'inf');
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
% (Proposition 20)
up_est = 2*exp(2*pi*c)*norm(K_2*K_1,'inf')/(exp(2*pi*c*N) - 1);

% Lower estimate for the LHS in Theorem 23
Nu = min(nu_(nu_>0));
low_est = rho_0^2*(1 + norm(K_3,'inf')*(Nu - C_2*pi/(2*m) - C_1)^(-1))^(-1);

% Warning if M,m or n are chosen too small and thus the lower estimate gets
% negative
if Nu < C_2*pi/(2*m)
    warning('m is chosen too small')
end
if Nu < C_1
    warning('M is chosen too small')
end