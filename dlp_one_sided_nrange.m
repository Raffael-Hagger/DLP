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
    M = 25;
end
if nargin > 5
    g = varargin{6};
else
    g = @(x) sin(pi*x).^2;
end
if nargin > 6
    s = varargin{7};
else
    s = 40;
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
B = Z + sum(exp(1i*J*t).*Pos + exp(-1i*J*t).*Neg,3);

% Computing T^{p,N}
e = exp(2i*pi*x');
T = zeros(2*p+1);
for j = 1:2*p+1
    for k = 1:2*p+1
        T(j,k) = (e').^(j-p-1)*B*e.^(k-p-1)/(N^2);
    end
end

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