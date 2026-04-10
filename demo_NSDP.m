%% Demo: Convex NSDP with ell_1-Regularization
% Problem formulation:
% \min \psi(x) = \sum_{i=1}^{n}( \frac{1}{4}d_i x_i^4 + \frac{1}{3}c_i |x_i|^3)
%                + \frac{1}{2}x^\top Q x + b^\top x + \rho \sum_{i=1}^n|x_i|
% s.t. G(x) = - A_0 - \sum_{i=1}^n x_i A_i \in \mathcal{S}_{-}^{m}.
% This demo studies the effect of different choices of $\bar r$ and $\bar s$, 
% which govern the rate of decrease of $\{\mu_{k}\}$

clear; clc
prob = 'a'; % problem

%% switch problem
switch prob
    case 'a'
        rand('seed', 2024)
        randn('seed', 2024)
        n = 500; 
        m = 500;
        rho = 1;

        % parameters
        rbars = [.33 0.6 .9]; 
        sbars = [0, 3];
        [obj, const, x0] = NSDP(m, n);

    case 'b'
        rand('seed', 2025)
        randn('seed', 2025)
        n = 1000; 
        m = 500;
        rho = 1;

        % parameters
        rbars = [.33 0.6 .9]; 
        sbars = [0, 3];
        [obj, const, x0] = NSDP(m, n);
end

%% print
% recorder
timestamp = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
fID = fopen('log.txt', 'a');
fprintf(fID, '  Ell-one regularized NSDP with (n,m) = (%3d,%3d), runing at %10s \n', n, m, char(timestamp));
fprintf(fID, '-------------------------------------------------------------------------------------------\n');
fprintf(fID, '%5s / %5s %10s   %9s   %8s   %5s   %7s   %4s   %8s   %8s\n', ...
             'rbar', 'sbar', 'f', 'g', 'time', 'iter', 'mean ls', 'mu0', 'mu', 'L_mu');
fprintf(fID, '-------------------------------------------------------------------------------------------\n');

fID_table = fopen('./tables.txt', 'a');
fprintf(fID_table, '  Ell-one regularized NSDP with (n,m) = (%3d,%3d), runing at %10s \n', n, m, char(timestamp));
fprintf(fID_table, '-------------------------------------------------------------------------------------------\n');
style_tex = {'{\color{blue}\bf - - - }', '{\color{blue}\bf $\cdots$ }', '{\color{blue}\bf --- }',...
    '{\color{red}\bf - - - }', '{\color{red}\bf $\cdots$ }', '{\color{red}\bf --- }'};
fprintf(fID_table, '                       tabale (%s)      \n', prob);
fprintf(fID_table, '----------------------------------------------------------------------\n');

%% sMBA loop
t = 0;
for ii = 1:length(sbars)
    sbar = sbars(ii);
for jj = 1:length(rbars)
    t = t+1;
    rbar = rbars(jj);
    Mths{t} = sprintf('$(%3.2g, %1d)$', rbar, sbar);
    fprintf('(rbar, sbar) = (%3.2g, %1d)...\n', rbar, sbar)

    % options
    opt = [];           
    opt.x0 = x0;
    opt.rbar = rbar; 
    opt.sbar = sbar;
    opt.epsilon = 1e-7;
    opt.isdisp = 1;
    opt.K = 5000;   
    opt.maxiter = 5000;

    % run sMBA
    tic_temp = tic;
    [xtemp, psi_temp, g_temp, out_temp] = ...
        sMBA_NSDP(obj, const, rho, opt);
    time_temp = toc(tic_temp);
    
    xs{t} = xtemp;
    psi_vals{t} = psi_temp;
    gvals{t} = g_temp;
    outs{t} = out_temp;
    times(t) = time_temp;

    % print results
    fprintf('%4s / %4s   %10s     %9s    %8s     %5s   %7s   %4s   %8s  %8s\n', ...
             'rbar', 'sbar', 'obj', 'g', 'time', 'iter', 'mean ls', 'mu0', 'mu', 'L_mu');
    fprintf('%4.2g / %1d | %10.10f | %9.2e | %8.5f | %5d | %7.2f | %4.2g | %8.2e | %8.2e\n', ...
        opt.rbar, opt.sbar, psi_temp, g_temp, time_temp, ...
        out_temp.k, mean(out_temp.ls(:,2),'all'), out_temp.mu0, out_temp.mu, out_temp.L_mu);
    fprintf('-------------------------------------------------------------------------------------------\n\n');
    
    fprintf(fID, '%4.2g / %1d | %10.10f | %9.2e | %8.5f | %5d | %7.2f | %4.2g | %8.2e | %8.2e\n', ...
        opt.rbar, opt.sbar, psi_temp, g_temp, time_temp, ...
        out_temp.k, mean(out_temp.ls(:,2),'all'), out_temp.mu0, out_temp.mu, out_temp.L_mu);

    fprintf(fID_table, '%27s & $(%4.2g, %1d)$ & $%12.5f$ & $%7.1f$\\\\\n',...
        style_tex{t}, rbars(jj), sbars(ii), psi_vals{t}, times(t));
end
end
fprintf(fID, '-------------------------------------------------------------------------------------------\n');

%% cvx solver
fprintf('cvx...\n')
st = tic;
cvx_solver SDPT3
cvx_begin SDP
    variable x_cvx(n)
    expression M(m,m)
    M = -const.A0 - reshape(x_cvx'*const.A, m, m);
    minimize( obj.d'*x_cvx.^4 / 4 + obj.c'*pow_abs(x_cvx, 3) / 3 + .5*quad_form(x_cvx, obj.Q) + obj.b' * x_cvx  + rho*norm(x_cvx, 1)) % 
    subject to
    M<=0;
cvx_end
time_cvx = toc(st);

% calculate lambda_{max}(G) 
G = -const.A0 - reshape(x_cvx'*const.A, m, m); % m x m
[U, S] = eig(G);
y = diag(S);
g_cvx = max(y);

% print results
disp('Method         obj              g        time');
fprintf('%6s | %10.8f | %2.2e | %8.5f\n', ...
    'cvx', cvx_optval, g_cvx, time_cvx);
fprintf('-------------------------------------------------------------------------------------------\n');

fprintf(fID, '%4s %3s | %10.10f | %9.2e | %8.5f\n', ...
    '', 'cvx', cvx_optval, g_cvx, time_cvx);
fprintf(fID, '-------------------------------------------------------------------------------------------\n\n');

fprintf(fID_table, '%28s &  %4s  & $%12.5f$ & $%7.1f$ \n', '', 'cvx', cvx_optval, time_cvx);
fprintf(fID_table, '----------------------------------------------------------------------\n\n');
fclose(fID);
fclose(fID_table);

%% figure
close all
figure; hold on
pos = [1 1 5 4.2];
set(gcf,'Units','Inches');
set(gcf,'Position',pos);
style = {'--b', ':b', '-b', '--r', ':r', '-r'};
tt = 6;

% plot sMBA
for t = 1:tt
    temp = (outs{t}.psi_vals-cvx_optval)/max(1,abs(cvx_optval));
    plot(log10(temp(2:end-1)), style{t});
end
axis([0, 5000, -9, 4.5])
title(sprintf('(%s) $(n, m) = (%d, %d)$', prob, n, m), 'Interpreter','latex')
xlabel('iteration $k$','Interpreter','latex')
ytickformat('10^{%.2g}')
ylabel('$\omega_{k} $', 'Interpreter','latex')


