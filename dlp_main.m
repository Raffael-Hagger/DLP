function[] = dlp_main

% Choosing different options to compute
estimates = 1;
spectrum = 1;
nrange = 1;
sides = 1; % set to 1 for the one-sided case, 2 for two-sided case

% Fixing the parameters, functions and norms
if sides == 1
    % General parameters
    a = 1/2;
    m = 400;
    M = 25;

    % Only needed for estimates
    c = 0.03;
    rho_0 = 1/2;

    % Only needed for spectrum and numerical range
    N = 256;

    % Only needed for numerical range
    t = 0;
    p = 10;
    s = 40;

    % Function g and its norm estimates
    g = @(x) sin(pi*x).^2; 
    ng = 1; % sup norm of g
    ng1 = pi; % sup norm of g'
    ng2 = 2*pi^2; % sup norm of g''
    nimg = (1/2)*sinh(2*pi*c); % c norm of Im(g)
    nimg1 = pi*sinh(2*pi*c); % c norm of Im(g')
    ngc = cosh(pi*c)^2; % c norm of g
    ngc1 = pi*cosh(2*pi*c); % c norm of g'
    ngc2 = 2*pi^2*cosh(2*pi*c); % c norm of g''
elseif sides == 2
    % General parameters
    a = 3/4;
    m = 400;
    M = 50;

    % Only needed for estimates
    c = 0.05;
    rho_0 = 1/2;
    
    % Only needed for spectrum and numerical range
    N = 32;

    % Only needed for numerical range
    t = 0;
    p = 10;
    s = 40;
    
    % Function g_+,g_- and their norm estimates
    g_p = @(x) 1;
    g_m = @(x) 1;
    ng_p = 1; % sup norm of g_+
    ng_m = 1; % sup norm of g_-
    ng1_p = 0; % sup norm of g_+'
    ng1_m = 0; % sup norm of g_-'
    ng2_p = 0; % sup norm of g_+''
    ng2_m = 0; % sup norm of g_+''
    nimg_p = 0; % c norm of Im(g_+)
    nimg_m = 0; % c norm of Im(g_-)
    nimg1_p = 0; % c norm of Im(g_+')
    nimg1_m = 0; % c norm of Im(g_-')
    ngc_p = 1; % c norm of g_+
    ngc_m = 1; % c norm of g_-
    ngc1_p = 0; % c norm of g_+'
    ngc1_m = 0; % c norm of g_-'
    ngc2_p = 0; % c norm of g_+''
    ngc2_m = 0; % c norm of g_-''
end

% If estimates is set to 1, this computes upper and lower estimates as well
% as the maximum spectral radius for N = 2^j
up_est = 1;
low_est = 0;
if estimates == 1
    j = 0;
    while up_est >= low_est
        j = j+1;
        if sides == 1
            [up_est,low_est,Rho] = dlp_test_one_sided(2^j,a,c,m,M,rho_0,g,ng,ng1,ng2,nimg,nimg1,ngc,ngc1,ngc2);
        elseif sides == 2
            [up_est,low_est,Rho] = dlp_test_two_sided(2^j,a,c,m,M,rho_0,g_p,g_m,ng_p,ng_m,ng1_p,ng1_m,ng2_p,ng2_m,nimg_p,nimg_m,nimg1_p,nimg1_m,ngc_p,ngc_m,ngc1_p,ngc1_m,ngc2_p,ngc2_m);
        end
        x(j) = 2^j;
        y1(j) = up_est;
        y2(j) = low_est;
        rho(j) = Rho;
    end

% Plot
    figure
    yyaxis left
    loglog(x,y1,'-sk',x,y2,'-^k')
    yyaxis right
    semilogx(x,rho,'-or')
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'r';
    xlabel('N')
    xticks(x)
end

% If spectrum is set to 1, this plots the spectrum for the chosen N
if spectrum == 1
    if sides == 1
        dlp_one_sided_spectrum(N,a,m,M,g)
    elseif sides == 2
        dlp_two_sided_spectrum(N,a,m,M,g_p,g_m)
    end
end

% If nrange is set to 1, this plots the (support lines of the) numerical
% range for the chosen t, p and N. s denotes the number of equally spaced
% support lines
if nrange == 1
    if sides == 1
        dlp_one_sided_nrange(N,a,t,p,M,g,s)
    elseif sides == 2
        dlp_two_sided_nrange(N,a,t,p,M,g_p,g_m,s)
    end
end
