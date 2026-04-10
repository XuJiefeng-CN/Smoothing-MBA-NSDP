function [obj, const, x0] = NSDP(m, n, lambda1, sp)
% This function generates the data of the following nonlinear SDP: 
% \min \psi(x) = \sum_{i=1}^{n}( d_i x_i^4/4 + c_i |x_i|^3/3) + x^\top Q x / 2 
%                   + b^\top x + rho \sum_{i=1}^n|x_i|
% s.t. G(x) = -A_0 - \sum_{i=1}^n x_i A_i \in \mathcal{S}_{-}^{m}.

if ~exist('lambda1', 'var')
    lambda1 = 100;
end
if ~exist('sp', 'var')
    sp = .2;   % Sparsity level
end

% Construct objective
[U, ~] = qr(randn(n));        % orthogonal matrix 
Lambda = lambda1*full(sprand(n,1, sp));
Lambda_sqrt = sqrt(Lambda);
Utemp = U*diag(Lambda_sqrt);

obj.Q = Utemp*Utemp'; % n x n
id = Lambda>0;
obj.b = 10+randn(n, 1);
obj.b = U*(id.*obj.b);
obj.c = lambda1*full(sprand(n,1, sp));
obj.d = lambda1*full(sprand(n,1, sp));

% Construct constraints
const.A =zeros(n, m, m);
for i=1:n
    [U, ~] = qr(randn(m));        % orthogonal matrix
    Lambda = lambda1*full(sprand(m, 1, sp));
    Lambda_sqrt = sqrt(Lambda);
    Utemp = U*diag(Lambda_sqrt);
    const.A(i, :, :) = Utemp*Utemp'; % m x m
end
const.A = reshape(const.A, n, []);
[U, ~] = qr(randn(m));        % orthogonal matrix
Lambda = lambda1*(.1 + .9*rand(m,1));
Lambda_sqrt = sqrt(Lambda);
Utemp = U*diag(Lambda_sqrt);
const.A0 = Utemp*Utemp'; % m x m
x0 = zeros(n, 1);
end








