function[L_c,R_c,Rho] = dlp_test_two_sided(varargin)

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
t = ((1:m) - 1/2)*pi/m;

% The estimates for \|\tilde{K}^{\pm}_t\|_{c,0},\|\tilde{L}^{\pm}_t\|_{c,0}
% and \|\tilde{K}^{\pm}_t\|_{0,c}, \|\tilde{L}^{\pm}_t\|_{0,c},
% respectively (Prop. 4.17)
I_c_p = nimg_p + nimg1_p/lga;
I_c_m = nimg_m + nimg1_m/lga;
k_0 = max(ng_p,ng_m)*lga/a;
k_c_p = max(ngc_p,ng_m)*lga/a;
k_c_m = max(ngc_m,ng_p)*lga/a;
F_0_p = ng_p + ng1_p/lga;
F_0_m = ng_m + ng1_m/lga;
F_c_p = ngc_p + ngc1_p/lga;
F_c_m = ngc_m + ngc1_m/lga;
G_0_p = ng1_p + ng2_p/lga;
G_0_m = ng1_m + ng2_m/lga;
G_c_p = ngc1_p + ngc2_p/lga;
G_c_m = ngc1_m + ngc2_m/lga;
if c > acos(a)/lga || I_c_p >= 1 || I_c_m >= 1 || nimg_p + a*c*k_c_p >= a^2 || nimg_m + a*c*k_c_m >= a^2
    warning('Conditions on c violated.')
end
K_p = (1 + F_c_p^2)^(1/4)/(pi*(1 - I_c_p^2))*(G_c_p*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + lga*(F_c_p + F_0_p)*BStar10(a));
K_m = (1 + F_c_m^2)^(1/4)/(pi*(1 - I_c_m^2))*(G_c_m*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + lga*(F_c_m + F_0_m)*BStar10(a));
L_p = (2*pi)^(-1)*(lga*F_0_m + k_c_p)/(1 - a^(-4)*(nimg_p + a*c*k_c_p)^2)*(1 + F_c_p^2)^(1/4)*CStar10(a);
L_m = (2*pi)^(-1)*(lga*F_0_p + k_c_m)/(1 - a^(-4)*(nimg_m + a*c*k_c_m)^2)*(1 + F_c_m^2)^(1/4)*CStar10(a);
C_5 = [K_m L_m;L_p K_p];
K_p = (1 + F_0_p^2)^(1/4)/(pi*(1 - I_c_p^2)^(5/4))*(G_c_p*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + 2*lga*F_c_p*BStar10(a));
K_m = (1 + F_0_m^2)^(1/4)/(pi*(1 - I_c_m^2)^(5/4))*(G_c_m*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + 2*lga*F_c_m*BStar10(a));
L_p = (2*pi)^(-1)*(lga*F_c_m + k_c_m)/(1 - a^(-4)*(nimg_m + a*c*k_c_m)^2)*(1 + F_0_p^2)^(1/4)/(1 - I_c_m^2)^(1/4)*CStar10(a);
L_m = (2*pi)^(-1)*(lga*F_c_p + k_c_p)/(1 - a^(-4)*(nimg_p + a*c*k_c_p)^2)*(1 + F_0_m^2)^(1/4)/(1 - I_c_p^2)^(1/4)*CStar10(a);
C_6 = [K_m L_m;L_p K_p];
C_4 = 2*exp(2*pi*c)*norm(C_6*C_5,'inf'); % (4.70)

% Estimates for \|\tilde{K}^{\pm}_t\|_{\infty} and
% \|\tilde{L}^{\pm}_t\|_{\infty} (see (4.69))
K_p = (1 + F_0_p^2)^(1/4)/pi*(G_0_p*(1 + a^(1/2) + a^(-1/2)/(4*a^2)) + 2*lga*F_0_p*BStar10(a));
K_m = (1 + F_0_m^2)^(1/4)/pi*(G_0_m*(1 + a^(1/2) + a^(-1/2)/(4*a^2)) + 2*lga*F_0_m*BStar10(a));
L_p = (2*pi)^(-1)*(lga*F_0_m + k_0)*(1 + F_0_p^2)^(1/4)*CStar10(a);
L_m = (2*pi)^(-1)*(lga*F_0_p + k_0)*(1 + F_0_m^2)^(1/4)*CStar10(a);
C_3 = norm([K_m L_m;L_p K_p],'inf'); % (4.69)

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

% Estimates for the constants C_1(M) and C_2 := \|B_N^M\|_{\infty} (see
% (4.66) and (4.68))
C_1_p = (1 + F_0_p^2)^(1/4)*(2*lga*F_0_p*log(tanh(abs((M-1)*log(a))/4))/(sqrt(a)*log(a)) + (lga*F_0_m + k_0)*2*atan(a^(M/2))/(a^2*abs(log(a))));
C_1_m = (1 + F_0_m^2)^(1/4)*(2*lga*F_0_m*log(tanh(abs((M-1)*log(a))/4))/(sqrt(a)*log(a)) + (lga*F_0_p + k_0)*2*atan(a^(M/2))/(a^2*abs(log(a))));
C_1 = max(C_1_p,C_1_m)/(2*pi);
C_2 = norm(sum([J J;J J].*abs(Pos) + [J J;J J].*abs(Neg),3),'inf');

% Summing the terms up
rho = zeros(m,1);
nu = zeros(m,2);
nu_ = zeros(m-1,2);
progress = waitbar(0,'Progress');
for k = 1:m
    A = Z + sum(exp(1i*[J J;J J]*t(k)).*Pos + exp(-1i*[J J;J J]*t(k)).*Neg,3);

% Computing the spectral radii of B_{t_k,N}^M for k = 0,...,m
    rho(k) = abs(eigs(A,1));

% Computing the the smallest singular value of (B_{t_k,N}^M - \lambda_l I),
% where \lambda_l is chosen adaptively
    phi = 0;
    nu(k,1) = 1/norm(inv(A - rho_0*eye(2*N)),'inf');
    l = 1;
    while phi < 2*pi
        phi = phi + nu(k,l)/(2*rho_0);
        lambda = rho_0*exp(1i*phi);
        nu(k,l+1) = 1/norm(inv(A - lambda*eye(2*N)),'inf');
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

% Computation of L_c and R_c in (4.43)
L_c = C_4/(exp(2*pi*c*N) - 1);
R = min(nu_(nu_>0));
R_c = rho_0^2*(1 + C_3*(R - C_1 - C_2*pi/(2*m))^(-1))^(-1);

% Warning if M, m or n are chosen too small so that R_c is
% negative
if R_c < 0
    warning('R_c is negative (N = %d)',N)
end