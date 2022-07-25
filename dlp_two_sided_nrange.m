function[] = dlp_two_sided_nrange(varargin)

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
    t = varargin{3};
else
    t = 0;
end
if nargin > 3
    p = varargin{4};
else
    p = 10;
end
if nargin > 4
    M = varargin{5};
else
    M = 50;
end
if nargin > 5
    g_p = varargin{6};
    g_m = varargin{7};
else
    g_p = @(x) 1;
    g_m = @(x) 1;
end
if nargin > 7
    s = varargin{8};
else
    s = 40;
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

% The 0-th term in the series of \tilde{K}_t
x = 1/(2*N):1/N:1-1/(2*N);
y = x;
[Y,X] = meshgrid(y,x);
Z_K_p = kernel_K_p(X,Y);
Z_K_m = kernel_K_m(X,Y);
Z_L_p = kernel_L_p(X,Y);
Z_L_m = kernel_L_m(X,Y);
for q = 1:N
    Z_K_p(q,q) = kernel_sing_p(x(q));
    Z_K_m(q,q) = kernel_sing_m(x(q));
end

% The j-th term, j = -M,...,-1,1,...,M
j = 1:M;
[Y,X,J] = meshgrid(y,x,j);
Pos_K_p = kernel_K_p(X+J,Y);
Pos_K_m = kernel_K_m(X+J,Y);
Pos_L_p = kernel_L_p(X+J,Y);
Pos_L_m = kernel_L_m(X+J,Y);
Neg_K_p = kernel_K_p(X-J,Y);
Neg_K_m = kernel_K_m(X-J,Y);
Neg_L_p = kernel_L_p(X-J,Y);
Neg_L_m = kernel_L_m(X-J,Y);

% Summing up the kernels
B_11 = Z_K_m + sum(exp(1i*J*t).*Pos_K_m + exp(-1i*J*t).*Neg_K_m,3);
B_12 = Z_L_m + sum(exp(1i*J*t).*Pos_L_m + exp(-1i*J*t).*Neg_L_m,3);
B_21 = Z_L_p + sum(exp(1i*J*t).*Pos_L_p + exp(-1i*J*t).*Neg_L_p,3);
B_22 = Z_K_p + sum(exp(1i*J*t).*Pos_K_p + exp(-1i*J*t).*Neg_K_p,3);

% Computing T^{p,N}
e = exp(2i*pi*x');
T_11 = zeros(2*p+1);
T_12 = zeros(2*p+1);
T_21 = zeros(2*p+1);
T_22 = zeros(2*p+1);
for j = 1:2*p+1
    for k = 1:2*p+1
        T_11(j,k) = (e').^(j-p-1)*B_11*e.^(k-p-1)/(N^2);
        T_12(j,k) = (e').^(j-p-1)*B_12*e.^(k-p-1)/(N^2);
        T_21(j,k) = (e').^(j-p-1)*B_21*e.^(k-p-1)/(N^2);
        T_22(j,k) = (e').^(j-p-1)*B_22*e.^(k-p-1)/(N^2);
    end
end
T = [T_11 T_12;T_21 T_22];

% Preparing the plot
figure
hold on
axis equal

% Computing the support lines of the numerical range
absc = zeros(s,1);
for j = 1:s
    phi = 2*pi*j/s;
    absc(j) = eigs(exp(-1i*phi)*T+exp(1i*phi)*T',1,'largestreal')/2;
    line([cos(phi)*absc(j)+sin(phi),cos(phi)*absc(j)-sin(phi)],[sin(phi)*absc(j)-cos(phi),sin(phi)*absc(j)+cos(phi)],'Color','black')
end
max(abs(absc))