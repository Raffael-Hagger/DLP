function[] = dlp_one_sided_spectrum(varargin)

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
    m = varargin{3};
else
    m = 400;
end
if nargin > 3
    M = varargin{4};
else
    M = 25;
end
if nargin > 4
    g = varargin{5};
else
    g = @(x) sin(pi*x).^2;
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

% Summing the terms up and computing the eigenvalues
progress = waitbar(0,'Progress');
eigA = zeros(1,N*m);
for k = 1:m
    A = Z + sum(exp(1i*J*t(k)).*Pos + exp(-1i*J*t(k)).*Neg,3);
    eigA((k-1)*N+1:k*N) = eig(A);
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