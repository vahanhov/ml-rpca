function [ D, imDims, L0, S0 ] = loadSynthetic( dims, r, varargin )
%% Generates and return synthetically generated data for experiments
%
%   Author: Vahan Hovhannisyan, 2017.

if length(dims) ~= 3 && length(dims) ~= 2
    error('Input should be of Width * Height * Numimages');
end

disp(['Loading synthetic data with rank ' num2str(r)]);

if ~isempty(varargin)
    incoh = varargin{1};
else
    incoh = 1.25;
end

if numel(varargin) > 1
    sparsity = varargin{2};
else
    sparsity = 0.5;
end

if numel(varargin) > 1
    degree = varargin{2};
else
    degree = 100;
end

rng(100);
if length(dims) == 3
    m = prod(dims(1:2));
    n = dims(3);
    imDims = dims(1:2);
else
    m = prod(dims(1));
    n = dims(2);
    imDims = [];
end

if incoh == 1 %% Generate the low rank part
    M = randn(m, r); % see RPCA paper for the distribution that you have to sample
    N = randn(r, n);
    L0 = M * N / r;
else %% Generate Low Rank part with certain incoherence parameter
    Sig = diag(r);
    fact = (1-(1/(m*n)))/(r-1);
    vari = 0.05;
    for i = 1 : r-1
    cent = 1-(i-1)*fact;
    Sig(i, i) = cent*(1 + vari * cent * randn); %1
    end
    Sig(r, r) = 1;%/(m*n); % set condition number

    zU = zeros(m, r);
    zV = zeros(n, r);
    rp = randperm(m);
    z_u = floor(m / incoh^2);
    zU(rp(1:z_u), :) = randn(z_u, r);
    rp = randperm(n);
    z_v = floor(n / incoh^2);
    zV(rp(1:z_v), :) = randn(z_v, r);
    zU = zU./repmat(sqrt(sum(zU .* zU, 1)), m, 1);
    zV = zV./repmat(sqrt(sum(zV .* zV, 1)), n, 1); % set incoherence
    L0 = zU * Sig * zV';
end

%% Make the singular values of L quadratically decreasing
[UL, SL, VL] = lansvd(L0, r, 'L');
range = 1:r;
sL = SL(1,1) ./ range.^2;
L0 = UL * diag(sL) * VL';


%% Generate the sparse part
S0 = sign(randn(m, n));

%% Make the singular values of S quadratically decreasing
[US, SS, VS] = svd(S0, 'econ');
range = 1:size(SS, 1);
sS = SS(1,1) ./ range.^2;
S0 = US * diag(sS) * VS';

%% Make S0 sparse
inds = rand(m, n) < sparsity;
S0(inds) = 0;

% error('FIX THE SINGULAR VALUES OF SPARSE COMPONENT!');

% p = degree / sqrt(m*n); % generate sparse component
% S0 = rand(m, n);
% s_zero_idx = S0 > p;
% s_nonzero_idx = find(S0 <= p);
% S0(s_zero_idx) = 0;
% S0(s_nonzero_idx) = (r / sqrt(m * n)) *...
%     (rand(length(s_nonzero_idx), 1) / 2 + .5) .* sign(randn(length(s_nonzero_idx), 1));

%% Add the low rank and sparse parts
D = L0 + S0;

disp(['Generated a random low rank + sparse matrix of ' num2str(size(D)) ' dimensions.']);

end

