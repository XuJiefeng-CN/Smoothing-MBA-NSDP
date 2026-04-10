function [x, psi, g, out] = sMBA_NSDP(obj, const, rho, opt)
% The sMBA for solving the following NSDP with ell-one regularization:
% \min \psi(x) = \sum_{i=1}^{n}( \frac{1}{4}d_i x_i^4 + \frac{1}{3}c_i |x_i|^3)
%                + \frac{1}{2}x^\top Q x + b^\top x + \rho \sum_{i=1}^n|x_i|
% s.t. G(x) = - A_0 - \sum_{i=1}^n x_i A_i \in \mathcal{S}_{-}^{m},
% 'obj' structure with b, c, d, Q
% 'const' structure with Ai, i=0,1,...,n
% opt.maxiter: maximum iteration
% opt: parameters

% out.status: 
% case 1: res1<epsilon && res2<epsilon
% case 2: g(x0)>=0
% case 2: Fail to line-search
% case 3: Exceeds the maximum iteration

if ~isfield(opt, 'n0') 
    opt.n0 = 300;
end
if ~isfield(opt, 'isdisp')
    opt.isdisp = 0;
end
if opt.isdisp
    st = tic;
end
if ~isfield(opt, 'rbar')
    opt.rbar = .33;
end
if ~isfield(opt, 'sbar')
    opt.sbar = 6;
end
if ~isfield(opt, 'K')
    opt.K = 5000;
end
if ~isfield(opt, 'maxiter')
    opt.maxiter = 10000;
end
if ~isfield(opt, 'mu0')
    opt.mu0 = .9;
end


% initial point
x = opt.x0;           
epsilon = opt.epsilon;
m = length(const.A0);

% line search parameters
tau1 = .01;   
tau2 = .01; 
L_check = 1e-8;
L_hat = 1e8;
Lf = 1;
Lg = 1;

% smoothing parameters
mu0 = opt.mu0;
n0 = opt.n0;
nu0 = 1/(10*opt.n0 + 1);
r0 = .01;
s0 = 0;
sbar = opt.sbar;
rbar = opt.rbar;
alpha4 = 1e-5;
K = opt.K;

% recorders
ls = zeros(K+1,2);    % total iteration of line-search
psi_s = [];         % record the values of objective f

% require g(x0)<0
G = -const.A0 - reshape(x'*const.A, m, m); % m x m
G = .5*(G + G');
[U, S] = eig(G);
y = diag(S);
g = max(y);

if g>=0
    disp('g(x0)>=0')
    psi = [];
    out = [];
    out.status = 2;
    return
end



% find mu0>0 s.t. g_mu0(x0)<.1*g(x0)   
while 1
    u = exp((y-g)/mu0);
    u_sum = sum(u);
    g_mu = g + mu0*log(u_sum) + alpha4*mu0; % g_mu0(x0)

    % break condition 
    if g_mu<= .1*g
        break
    end

    % decrease mu0
    mu0 = .5*mu0;
end
mu = mu0;

% psi(x0)
Qx = obj.Q*x;
dx = obj.d.*x.^3;
cx = obj.c.*x.^2.*sign(x);
f = x'*(dx/4 + cx/3 + Qx/2 + obj.b);
df = dx + cx + Qx + obj.b;
psi = f + rho*sum(abs(x));
psi_s = [psi_s; psi];

%  nabla g_mu0(x0)
dh_mu = u/u_sum;            
dh_mu = U*(dh_mu.*U'); % nabla h_mu(G(x))
dg_mu = -const.A*reshape(dh_mu, [], 1);

if opt.isdisp
    % print
    fprintf('\n------------------------------- sMBA -----------------------------------------\n');
    fprintf('%5s   %5s   %5s   %10s   %9s   %8s   %8s   %8s   %8s   %8s    %8s   %8s   %8s   %8s   %8s\n', ...
            'k', 'is', 'js', 'psi', 'g', 'mu', 'Lf0', 'Lf1', 'Lg0', 'Lg1', 'time', 'res1', 'res2', 'lambda', 'Delta');
    fprintf('---------------------------------------------------------------------------------\n');
end
k = 0;
while k<=opt.maxiter
    % line-search
    i = 0;
    ls_suc = 0;
    if opt.isdisp
        Lf0 = Lf;
        Lg0 = Lg;
    end
    for j = 0:40
        Lmu = Lg/mu;
        % Solve subproblem
        RR = (dg_mu/Lmu)'*(dg_mu/Lmu) - 2*g_mu/Lmu;
        Delta = sqrt(dg_mu'*dg_mu - 2*Lmu*g_mu);
        xbar = x - 1/Lf * df;
        xhat = x - 1/Lmu * dg_mu;
        [x1, lambda1] = SubP_alpha(xbar, xhat, RR, Lf/rho);
        lambda1 = 2*lambda1/Lmu;

        % calculate Gx1 and gx1
        G = -const.A0 - reshape(x1'*const.A, m, m); % m x m
        G = .5*(G + G');
        [U, S] = eig(G);
        y = diag(S);
        g = max(y);
        
        % calculate g1_mu := g_mu(x1)
        u = exp((y-g)/mu);
        u_sum = sum(u);
        g1_mu = g + mu*log(u_sum) + alpha4*mu;
        
        if g1_mu<=0      % is feasible
            % calculate psi(x1)
            Qx = obj.Q*x1;
            dx = obj.d.*x1.^3;
            cx = obj.c.*x1.^2.*sign(x1);
            f1 = x1'*(dx/4 + cx/3 + Qx/2 + obj.b);
            psi1 = f1 + rho*sum(abs(x1));
            s = x1 - x;
            ss = s'*s;

            % Sufficient decreasing condition
            if psi1 <= psi - .5*(tau1 + tau2*lambda1/mu)*ss
                % Both conditions hold
                ls_suc = 1;
                break
            else
                % increase Lf
                Lf = 2*Lf;
                i = i+1;
            end
        end
        % increase Lg
        Lg = 2*Lg;
    end
    ls(k+1, 1) = i;
    ls(k+1, 2) = j;
    if opt.isdisp
        Lf1 = Lf;
        Lg1 = Lg;
    end

    % termination conditions
    res1 = sqrt(tau1*mu + tau2*lambda1)*sqrt(ss)/mu/max(1,norm(x1));
    res2 = - lambda1*(y'*u)/sum(u)/max(1,norm(x1));
    % dhx1mu0 = u/sum(u);
    % dhx1mu0 = U*(dhx1mu0.*U');
    % v = lambda1*dhx1mu0;
    % res3 = - (v(:)' * G(:))/max(1,norm(x1)); % slackness
    % if abs(res2 - res3)>1e-14
    %     abs(res2 - res3)
    % end

    if ~ls_suc || (res1<=epsilon && res2<=epsilon)
        % update x and psi
        if ls_suc
            x = x1;
            psi = psi1;
            psi_s = [psi_s; psi];
            out.status = 1;
        else
            out.status = 3;
            disp('Fail to line-search!!!')
        end
        if opt.isdisp
            ttime = toc(st);
            % Print
             fprintf('%5d | %5d | %5d | %10.9e | %9.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e |%8.3g | %4.2e | %4.2e | %4.2e | %4.2e\n', ...
                    k, i, j, psi, g, mu, Lf0, Lf1, Lg0, Lg1, ttime, res1, res2, lambda1, Delta);

        end
        break
    end

    % update mu(k+1)
    k1 = mod(k+1, n0+1);           % remainder
    kprod = k+1 - k1;
    kbar = kprod + nu0*k1;
    r = r0 + min(1, kbar/K)*(rbar-r0);     % r(k+1)
    s_hat = s0 + min(1, kbar/K)*(sbar - s0);
    mu = mu0*((kbar + 1).^(-r))./(log(3+kbar).^s_hat);

    % calculate dg_mu1
    u = exp((y-g)/mu);
    u_sum = sum(u);
    g_mu = g + mu*log(u_sum) + alpha4*mu;
    dh_x1mu1 = u/u_sum;
    dh_x1mu1 = U*(dh_x1mu1.*U'); % nabla h_mu
    dg_x1mu1 = -(const.A*reshape(dh_x1mu1, [], 1));
    delta_g = mu*(dg_x1mu1 - dg_mu);
    Lg = max(L_check, Lg/2);

    % calculate BB stepsize: Lg_bb2
    if sqrt(abs(delta_g'*s))>1e-12
        Lg_bb = (delta_g'*delta_g)/abs(delta_g'*s); % BB stepsize of gmu          
        if Lg_bb >= L_check && Lg_bb <= L_hat
            Lg = Lg_bb;
        end
    end

    df1 = dx + cx + Qx + obj.b;
    delta_f = df1 - df;
    Lf = max(L_check, Lf/2);

    % calculate BB stepsize: Lf_bb1
    if sqrt(ss)>1e-12
        Lf_BB = abs(delta_f'*s/ss);      % BB stepsize of f
        if Lf_BB >= L_check && Lf_BB <= L_hat
            Lf = Lf_BB;
        end
    end
    
    % calculate BB stepsize: Lf_bb2
    % if sqrt(abs(delta_f'*s))>1e-12
    %     Lf_BB = delta_f'*delta_f / abs(delta_f'*s);      % BB stepsize of f
    %     if Lf_BB >= L_check && Lf_BB <= L_hat
    %         Lf = Lf_BB;
    %     end
    % end    
    
    % update x, psi and g
    x = x1;
    psi = psi1;
    df = df1;
    dg_mu = dg_x1mu1;
    psi_s = [psi_s; psi];
    if opt.isdisp
        ttime = toc(st);
        % Print
        if mod(k, ceil(K/20))==0 || k<=10 || j>25
             fprintf('%5d | %5d | %5d | %10.9e | %9.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e |%8.3g | %4.2e | %4.2e | %4.2e | %4.2e\n', ...
                        k, i, j, psi, g, mu, Lf0, Lf1, Lg0, Lg1, ttime, res1, res2, lambda1, Delta);
        end
    end
    k = k + 1;
end
if k>opt.maxiter
    out.status = 4;
end
if opt.isdisp
    fprintf('---------------------------------------------------------------------------------\n\n');
end

% output
out.k = k;
out.mu0 = mu0;
out.mu = mu;
out.L_mu = Lmu;
out.ls = ls;
out.psi_vals = psi_s;
out.res_const = norm(max(G,0));
out.res1 = res1;
out.res2 = res2;
end