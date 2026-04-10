%% Demo: Convex NSDP with ell_1-Regularization
% Problem formulation:
% \min \psi(x) = \sum_{i=1}^{n}( \frac{1}{4}d_i x_i^4 + \frac{1}{3}c_i |x_i|^3)
%                + \frac{1}{2}x^\top Q x + b^\top x + \rho \sum_{i=1}^n|x_i|
% s.t. G(x) = - A_0 - \sum_{i=1}^n x_i A_i \in \mathcal{S}_{-}^{m}.
% This demo tests sMBA on a range of dimensions $m$ and $n$ with different
% $\epsilon$.
clear;

% large size
rand('seed', 2026)
randn('seed', 2026)
ms = [100, 500];
ns = [100, 500, 1000]; 
KK = 30; % repeated number

% small size
% ms = [20, 50];
% ns = [10, 50, 100]; 
% KK = 3; % repeated number



%% Set up results recorder
timestamp = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
fID = fopen('./log.txt', 'a');
fprintf(fID, '  Ell-one regularized NSDP runing at %10s \n',... 
    char(timestamp));
fprintf(fID, '-------------------------------------------------------------------------------------------\n');

tb_ID = fopen('./tables.txt', 'a');
fprintf(tb_ID, '  Ell-one regularized NSDP runing at %10s \n',... 
    char(timestamp));
fprintf(tb_ID, '-------------------------------------------------------------------------------------------\n');

%%
for ii = 1:length(ms)
m = ms(ii);
fprintf(tb_ID, '\\midrule\n');
fprintf(tb_ID, '\\multirow{%1d}{*}{\\tabincell{c}{ $%3d$    }}\n', length(ns), m);
for jj = 1:length(ns)
for kk = 1:KK
    t = 0;
    n = ns(jj);
    [obj, const, x0] = NSDP(m, n);
    rho = 1; % coefficient of ell-one norm regularization
    fprintf(fID, '%3d & %3d\n', m, n);
    fprintf('%3d & %3d\n', m, n);
    fprintf(fID, 'Method  %10s   %9s    %8s     %5s   %7s   %4s   %8s  %8s  %8s  %8s   %8s\n', ...
             'obj', 'g', 'time', 'iter', 'mean ls', 'mu0', 'mu', 'L_mu', 'res1', 'res2', 'status');
    fprintf('Method %10s   %9s    %8s     %5s   %7s   %4s   %8s  %8s  %8s  %8s\n', ...
             'obj', 'g', 'time', 'iter', 'mean ls', 'mu0', 'mu', 'L_mu', 'res1', 'res2');
    
    %% sMBA-1
    t = t+1;
    Methods{t} = 'sMBA-1';
    fprintf('sMBA...\n')
    % options
    opt = [];           
    opt.x0 = x0;
    opt.rbar = .9;
    opt.sbar = 3;
    opt.epsilon = 1e-5;
    opt.isdisp = 1;
    opt.K = 5000;   
    opt.maxiter = 5000;

    % run sMBA
    tic_temp = tic;
    [xtemp, psi_temp, g_temp, out_temp] = ...
        sMBA_NSDP(obj, const, rho, opt);
    time_temp = toc(tic_temp);
    
    %  save results
    xs{t, ii, jj, kk} = xtemp;
    psi_vals(t, ii, jj, kk) = psi_temp;
    gvals(t, ii, jj, kk) = g_temp;
    outs{t, ii, jj, kk} = out_temp;
    iters(t, ii, jj, kk) = out_temp.k;
    times(t, ii, jj, kk) = time_temp;

    fprintf(fID, '%5s | %10.10f | %9.2e | %8.5f | %5d | %7.2f | %4.2g | %8.2e | %8.2e | %8.2e | %8.2e | %1d\n', ...
        Methods{t}, psi_temp, g_temp, time_temp, ...
        out_temp.k, mean(out_temp.ls(:,2),'all'), out_temp.mu0, out_temp.mu, out_temp.L_mu, out_temp.res1, out_temp.res2, out_temp.status);
    
    fprintf('%5s | %10.10f | %9.2e | %8.5f | %5d | %7.2f | %4.2g | %8.2e | %8.2e | %8.2e | %8.2e | %1d\n', ...
        Methods{t}, psi_temp, g_temp, time_temp, ...
        out_temp.k, mean(out_temp.ls(:,2),'all'), out_temp.mu0, out_temp.mu, out_temp.L_mu, out_temp.res1, out_temp.res2, out_temp.status);
   
    %% sMBA-2
    t = t+1;
    Methods{t} = 'sMBA-2';
    fprintf('sMBA...\n')
    % options
    opt = [];           
    opt.x0 = x0;
    opt.rbar = .9;
    opt.sbar = 3;
    opt.epsilon = 1e-7;
    opt.isdisp = 1;
    opt.K = 5000;   
    % Maximum number of iterations
    opt.maxiter = 5000;

    % run sMBA
    tic_temp = tic;
    [xtemp, psi_temp, g_temp, out_temp] = ...
        sMBA_NSDP(obj, const, rho, opt);
    time_temp = toc(tic_temp);
    
    %  save results
    xs{t, ii, jj, kk} = xtemp;
    psi_vals(t, ii, jj, kk) = psi_temp;
    gvals(t, ii, jj, kk) = g_temp;
    outs{t, ii, jj, kk} = out_temp;
    iters(t, ii, jj, kk) = out_temp.k;
    times(t, ii, jj, kk) = time_temp;

    fprintf(fID, '%5s | %10.10f | %9.2e | %8.5f | %5d | %7.2f | %4.2g | %8.2e | %8.2e | %8.2e | %8.2e | %1d\n', ...
        Methods{t}, psi_temp, g_temp, time_temp, ...
        out_temp.k, mean(out_temp.ls(:,2),'all'), out_temp.mu0, out_temp.mu, out_temp.L_mu, out_temp.res1, out_temp.res2, out_temp.status);
    
    fprintf('%5s | %10.10f | %9.2e | %8.5f | %5d | %7.2f | %4.2g | %8.2e | %8.2e | %8.2e | %8.2e | %1d\n', ...
        Methods{t}, psi_temp, g_temp, time_temp, ...
        out_temp.k, mean(out_temp.ls(:,2),'all'), out_temp.mu0, out_temp.mu, out_temp.L_mu, out_temp.res1, out_temp.res2, out_temp.status);
    
     %% cvx solver
    t = t + 1;
    Methods{t} = 'CVX';
    fprintf('cvx...\n')
    st = tic;
    cvx_solver SDPT3
    cvx_begin SDP
        variable x_cvx(n)
        expression M(m,m)
        M = -const.A0 - reshape(x_cvx'*const.A, m, m);
        minimize( obj.d'*x_cvx.^4 / 4 + obj.c'*pow_abs(x_cvx, 3) / 3 + .5*quad_form(x_cvx, obj.Q) + obj.b' * x_cvx  + rho*norm(x_cvx, 1))
        subject to
        M<=0;
    cvx_end
    time_cvx = toc(st);

    % constraint violation
    G = -const.A0 - reshape(x_cvx'*const.A, m, m); % m x m
    G = .5*(G + G');
    [U, S] = eig(G);
    y = diag(S);
    g_cvx = max(y);

    % save results
    xs{t, ii, jj, kk} = x_cvx;
    psi_vals(t, ii, jj, kk) = cvx_optval;
    gvals(t, ii, jj, kk) = g_cvx;
    times(t, ii, jj, kk) = time_cvx;

    fprintf(fID, ' %5s | %10.10f | %2.2e | %8.5f\n', ...
        'cvx', cvx_optval, g_cvx, time_cvx);
    fprintf(fID, '-------------------------------------------------------------------------------------------\n\n');

    fprintf('  %5s | %10.10f | %2.2e | %8.5f\n', ...
        'cvx', cvx_optval, g_cvx, time_cvx);
    fprintf('-------------------------------------------------------------------------------------------\n\n');

end
% print table
tt = t;
fprintf('& %3d\n', n);
fprintf(tb_ID, '& %3d\n', n);
for t = 1:tt
    psi_temp = mean(psi_vals(t, ii, jj, :));
    g_temp = mean(max(gvals(t, ii, jj, :), 0));
    time_temp = mean(times(t, ii, jj, :));
    if t<tt
        iter_temp = mean(iters(t, ii, jj, :));
        out_temp = outs{t, ii, jj};
        fprintf('  & %5.4g / %8.8g / %2.2g / %5.1f\n', ...
            iter_temp, psi_temp, g_temp, time_temp);
        fprintf(tb_ID, '  & %5.4g / %8.8g / %2.2g / %5.1f\n', ...
            iter_temp, psi_temp, g_temp, time_temp);
    else
        fprintf('   &        %8.8g / %2.2g / %5.1f\\\\\n', psi_temp, g_temp, time_temp);
        fprintf(tb_ID, '   &        %8.8g / %2.2g / %5.1f\\\\\n', psi_temp, g_temp, time_temp);
    end
end
end
end
fclose(fID);
fclose(tb_ID);
