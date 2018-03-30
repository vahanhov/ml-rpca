function [ output ] = ML_FW_T(par)
%% This function is derived from FW_T  
%   by replacing the singular value thresholing part by 
%   a multilevel singular value thresholding
%
%   Author: Vahan Hovhannisyan, 2017.

disp([' == running with ' num2str(par.opts.depth) ' levels == ']);

starttime = par.starttime;

%% setup
M = par.D; % data matrix
[m,n] = size(M);
normM = norm(M, 'fro');
% if isfield(par, 'GT')
%     normGTL = norm(par.GT.L, 'fro');
%     normGTS = norm(par.GT.S, 'fro');
% end

% stopping criterion by default (if no specification)
epsilon = 10^-3; 
maxtime = par.maxtime;
multilevel = true;

if isfield(par, 'epsilon') epsilon = par.epsilon; end
if isfield(par, 'maxtime') maxtime = par.maxtime; end
if isfield(par, 'multilevel') multilevel = par.multilevel; end

lambda_1 = par.lambda_1;% * 0.9;
lambda_2 = par.lambda_2;% / 0.9;
iter = par.iter;
Omega = par.Omega; % obeational pattern

if isfield(par, 'opts')
    opts = par.opts; % Options for multilevel iterations
else
    opts.depth = 2;
end

flag_T = false;
if m > n
    
    M = M';
    Omega = Omega';
    [m,n] = size(M);
    flag_T = true;
    
    if isfield(par, 'GT')
        par.GT.L = par.GT.L';
        par.GT.S = par.GT.S';
    end
end

%% FW-T method

%% initialization
L = zeros(m,n); S = zeros(m,n);
t_1 = 0; t_2 = 0;

U_1 = 0.5*normM^2/lambda_1; % initial rough guess
U_2 = 0.5*normM^2/lambda_2; % initial rough guess

history = 0; % to store the values of each iteration
gt_error.L = zeros(iter, 1);
gt_error.S = zeros(iter, 1);
change.L = zeros(iter, 1);
change.S = zeros(iter, 1);


fprintf('%6s    %3s             %10s     %10s   %5s \n', ...
    'iter.', 'obj. ', 'rel. err.', 'decrease', 'count');

grad = L + S - M;
temp = Omega .* grad;  % gradient
history(1) = norm(temp,'fro')^2/2+lambda_1*t_1+lambda_2*t_2;
fprintf('|  %2d |   %10.5d   |   %7.5d       |  %7.5d   | %2.1d | \n', ...
            0,      history(1),     inf,   inf,      0);
count = 0;

%% loop
for k = 1 : iter
    
    %------------------linearization subproblem-----------------------%
    if multilevel
        [ D_L, U, ev, V ] = mlLinNuc( temp, [], opts );
        evs(k) = ev;
    else
        [U,ev,V] = lansvd(temp,1,'l'); % top singular pair    

        D_L = -U*V';
    end
    
    if lambda_1 >= ev
        V_L = zeros(m,n); V_t_1 = 0;
    else
        V_L = U_1*D_L; V_t_1 = U_1;
    end
    
    % [mag ind] = max(reshape(abs(temp), numel(abs(temp)), 1));
    [mag, ind] = max(abs(temp(:)));
    j = floor((ind-1)/m)+1; i = mod(ind-1,m)+1;
    sign_ = sign(temp(i,j));
    D_S = zeros(m,n);
    D_S(i,j) = -sign_;
    
    if lambda_2 >= mag
        V_S = 0; V_t_2 = 0;
    else
        V_S = U_2*D_S; V_t_2 = U_2;
    end
    
    
    %--------------------- use QP (exact search) ---------------------%
    H = zeros(2,2);
    temp_1 = Omega.*(V_L-L); temp_2 = Omega.*(V_S-S); % temp = L + S - M;
    H(1,1) = norm(temp_1, 'fro')^2;
    H(2,2) = norm(temp_2, 'fro')^2;
    H(1,2) = temp_1(:)'*temp_2(:);
    H(2,1) = H(1,2);
    f = zeros(1,2);
    f(2) = temp(:)'*temp_2(:);
    f = f + [lambda_1*(V_t_1-t_1), lambda_2*(V_t_2-t_2)];
    
    % using QP solvers
    lb = zeros(2,1);
    ub = [1;1];
    f(1) = temp(:)'*temp_1(:);
    options.Display = 'off';
    options.TolFun = 1e-12;
    x = quadprog(H,f,[],[],[],[],lb,ub,[],options);    
    
    %----------------------- update L and S --------------------------%
    L_old = L;
    S_old = S;
    alpha = x(1);
    beta = x(2);
    L = (1-alpha)*L + alpha*V_L; t_1 = t_1 + alpha*(V_t_1-t_1);
    S = (1-beta)*S + beta*V_S; t_2 = t_2 + beta*(V_t_2-t_2);
    
    change.L(k) = norm(L - L_old, 'fro');
    change.S(k) = norm(S - S_old, 'fro');
    
    % Log the error from ground truth
    if isfield(par, 'GT')
        if isfield(par.GT, 'L') && ~isempty(par.GT.L)
            gt_error.L(k) = norm(par.GT.L - L, 'fro');
        end
        if isfield(par.GT, 'S') && ~isempty(par.GT.S)
            gt_error.S(k) = norm(par.GT.S - S, 'fro');
        end
    end
    
    %------------------------ thresholding----------------------------%
    temp_3 = S-Omega.*(L+S-M);
    S = max(temp_3 - lambda_2, 0);
    S = S + min(temp_3 + lambda_2, 0);
    t_2 = norm(S(:),1);
    
    %----------- update U_1 and U_2 to a better esitmate -------------%
	grad = L + S - M;
    temp = Omega .* grad;  % gradient
    history(k+1) = norm(temp,'fro')^2/2+lambda_1*t_1+lambda_2*t_2;
    U_1 = min(history(k+1)/lambda_1, U_1);
    U_2 = min(history(k+1)/lambda_2, U_2);
    
    decrease = history(k)-history(k+1);
    rel_err =  abs(decrease)/history(k);
    
    if rel_err <epsilon
        count = count+1;
    else
        count=0;
    end
    if count == 5
        fprintf('|  %2d |   %10.5d   |   %10.5d   |%10.5d | %2.1d | \n', ...
            k, history(k+1), rel_err, decrease, count);
        break;
    end
    
    if par.logs
        fprintf('|  %2d |   %10.5d   |   %10.5d   |%10.5d | %2.1d | \n', ...
                    k, history(k+1), rel_err, decrease, count);
    end
    
    if maxtime > 0 && (cputime - starttime > maxtime)
        disp(['=== Maximum time reached (' num2str(maxtime) ' seconds)']);
        break
    end
end

if flag_T  
   
    % transpose back
    M = M';
    L = L';
    S = S';
    Omega = Omega';
    
end


%% output
output.L = L;
output.S = S;
output.hist = history;
output.iter = k;

output.gt_error.L = gt_error.L(1:k);
output.gt_error.S = gt_error.S(1:k);

output.change.L = change.L(1:k);
output.change.S = change.S(1:k);


end