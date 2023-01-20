function[] = dlp_one_sided_nrange(varargin)

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
    t = varargin{4};
else
    t = 0;
end
if nargin > 4
    p = varargin{5};
else
    p = 10;
end
if nargin > 5
    M = varargin{6};
else
    M = 25;
end
if nargin > 6
    n = varargin{7};
else
    n = 40;
end
if nargin > 7
    g = varargin{8};
    ng = varargin{9}; % sup norm of g
    ng1 = varargin{10}; % sup norm of g'
    ng2 = varargin{11}; % sup norm of g''
    nimg = varargin{12}; % c norm of Im(g)
    nimg1 = varargin{13}; % c norm of Im(g')
    ngc = varargin{14}; % c norm of g
    ngc1 = varargin{15}; % c norm of g'
    ngc2 = varargin{16}; % c norm of g''
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

% The 0-th term in the series of \tilde{K}_t
x = 1/(2*N):1/N:1-1/(2*N);
y = x;
[Y,X] = meshgrid(y,x);
Z = kernel(X,Y);
for q = 1:N
    Z(q,q) = kernel_sing(x(q));
end

% The j-th term, j = -M,...,-1,1,...,M
j = 1:M;
[Y,X,J] = meshgrid(y,x,j);
Pos = kernel(X+J,Y);
Neg = kernel(X-J,Y);

% Summing up the kernel
A = Z + sum(exp(1i*J*t).*Pos + exp(-1i*J*t).*Neg,3);

% Computing T^{p,t,N,M} (see 4.74)
e = exp(2i*pi*x');
T = zeros(2*p+1);
for j = 1:2*p+1
    for k = 1:2*p+1
        T(j,k) = (e').^(j-p-1)*A*e.^(k-p-1)/(N^2);
    end
end

% Computing the largest eigenvector \lambda_l of
% \Re(e^{-i\theta_l}T^{p,t,N,M}) and a corresponding eigenvector x_l for
% l = 1,...,n in order to get the points z_l
z = zeros(1,n+1);
for l = 0:n-1
    theta = 2*pi*l/n;
    [x,~] = eigs(exp(-1i*theta)*T+exp(1i*theta)*T',1,'largestreal');
    z(l+1) = x'*T*x;
end
z(n+1) = z(1);

% The estimates for C_1, \|\tilde{K}_t\|_{c,0}, \|\tilde{K}_t\|_{0,c} and
% C_7 (see (4.34), Prop. 4.8 and Cor. 4.24)
I_c = nimg + nimg1/lga;
F_0 = ng + ng1/lga;
F_c = ngc + ngc1/lga;
G_c = ngc1 + ngc2/lga;
if c > acos(a)/lga || I_c >= 1
    warning('Condition (4.21) violated.')
end
C_1 = (2*lga*F_0/pi)*(1 + F_0^2)^(1/4)*log(tanh(abs((M-1)*log(a))/4))/(sqrt(a)*log(a));
C_5 = (1 + F_c^2)^(1/4)/(pi*(1 - I_c^2))*(G_c*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + lga*(F_c+F_0)*BStar10(a));
C_6 = (1 + F_0^2)^(1/4)/(pi*(1 - I_c^2)^(5/4))*(G_c*(1 + a^(1/2) + a^(-1/2))/(4*a^2) + 2*lga*F_c*BStar10(a));
C_7 = 2*(2*p+1)*exp(pi*(2*p+1)*c)*(C_5+C_6)/(exp(2*pi*c*N) - 1) + C_1;

% A lower estimate for the numerical range (see Cor. 4.25)
theta_max = max(diff(angle(z)));
R_min = min(abs(z));
RStar = R_min*cos(theta_max/2) - C_7;

% Preparing the plot
% figure
hold on
axis equal
plot(real(z),imag(z),'-k')
plot(RStar*exp(2i*pi*(0:0.001:1)),'-b')
