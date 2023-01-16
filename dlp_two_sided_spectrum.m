function[] = dlp_two_sided_spectrum(varargin)

% Setting some default options
if nargin > 0
    N = varargin{1};
else
    N = 16;
end
if nargin > 1
    a = varargin{2}; % a := \alpha
else
    a = 7/8;
end
if nargin > 2
    m = varargin{3};
else
    m = 2000;
end
if nargin > 3
    M = varargin{4};
else
    M = 200;
end
if nargin > 4
    g_p = varargin{5};
    g_m = varargin{6};
else
    g_p = @(x) 1;
    g_m = @(x) 1;
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

% Summing the terms up and computing the eigenvalues
progress = waitbar(0,'Progress');
eigA = zeros(1,2*N*m);
for k = 1:m
    A = Z + sum(exp(1i*[J J;J J]*t(k)).*Pos + exp(-1i*[J J;J J]*t(k)).*Neg,3);
    eigA((k-1)*2*N+1:k*2*N) = eig(A);
    text = strcat('Progress:',32,'N',32,'=',32,sprintf('%1.0f',N),44,32,'k',32,'=',32,sprintf('%1.0f',k));
    waitbar(k/m,progress,text)
end
delete(progress);

% Preparing the plot
figure
hold on
axis equal
plot(real(eigA),imag(eigA),'.k')
plot(real(eigA),-imag(eigA),'.k')
plot(0,0,'.k')