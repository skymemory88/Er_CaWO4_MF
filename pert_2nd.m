function results = pert_2nd(H0, H_pert, varargin)
% pert_2nd() - Second-order perturbation theory for arbitrary Hamiltonians
%
% Calculates first- and second-order perturbation corrections for a quantum
% system with an arbitrary unperturbed Hamiltonian H0 and one or more
% perturbation operators {V_i}. The routine handles non-degenerate spectra
% and (by default) performs proper *degenerate* PT within each degenerate
% manifold for first order, and excludes intra-manifold couplings in the
% second-order sums. It can also compute the full second-order tensor K_ij
% when the spectrum is non-degenerate. For degenerate spectra, supply a
% specific direction (linear combination) of perturbations via 'direction'
% to obtain unambiguous per-state second-order coefficients.
%
% SYNTAX:
%   results = perturb2nd(H0, V)                            % single perturbation
%   results = perturb2nd(H0, {Vx,Vy,Vz})                   % vector perturbations
%   results = perturb2nd(H0, H_pert, 'direction', uhat)    % along direction uhat
%   results = perturb2nd(H0, H_pert, 'calcTensor', true)   % full K_ij (non-deg only)
%
% INPUTS:
%   H0      - [N x N] Hermitian unperturbed Hamiltonian
%   H_pert  - Either a single [N x N] perturbation matrix V, or a cell array
%             {V1, V2, ... , Vk} of perturbation components
%
% OPTIONS:
%   'threshold'    - Absolute degeneracy threshold for energies (default: 1e-12)
%                    Levels |E_n - E_m| < threshold are treated as degenerate
%   'verbose'      - Logical, print notes (default: false)
%   'calcTensor'   - Logical, compute full second-order tensor K_ij per level
%                    (size [N x k x k]). Only meaningful if the spectrum is
%                    non-degenerate. If degeneracies are detected, K_ij is
%                    withheld and an effective block tensor is returned instead.
%                    Default: true
%   'calcWaveFunc' - Logical, compute first-order wavefunction corrections
%                    (default: true). Returned as cell array, one [N x N] matrix
%                    of coefficients per perturbation component
%   'direction'    - Real vector u of length k. If provided, the code computes
%                    the *combined* perturbation V_u = sum_i u_i V_i and returns
%                    first- and second-order corrections for this single effective
%                    perturbation (no K_ij). This is the recommended way to get
%                    unambiguous results in the presence of degeneracies (e.g. at
%                    B=0 for Kramers doublets)
%
% OUTPUTS:
%   results - Structure containing:
%     .eigE0        - [N x 1] Sorted unperturbed energies
%     .eigW0        - [N x N] Corresponding eigenvectors (columns)
%     .Hproj        - {k} Projected perturbations W0' * V_i * W0 (Hermitianized)
%     .dE1          - [N x k] First-order energy corrections (see notes)
%     .dE2          - [N x k] Second-order diagonal coefficients (if calcTensor=false)
%                     OR [N x k x k] K_ij tensor (if calcTensor=true and non-deg)
%     .dPsi1        - {k} First-order wavefunction corrections (if requested)
%     .hasDegeneracy- logical True if any degenerate manifold found under 'threshold'
%     .blocks       - cell Cell array of index vectors for each degenerate block
%     .blockK       - struct (Only if degeneracy && calcTensor) block-level second-order
%                     tensors K_ij for each degenerate manifold
%     .parameters   - struct Echo of options used
%
% THEORY REMINDERS:
%   E_n^(1)(i)   = <n|V_i|n>  (Hellmann–Feynman; for degeneracies, diagonalize V_i in block)
%   E_n^(2)(i)   = sum_{m!=n} |V_i(nm)|^2 / (E_n - E_m)
%   K_ij^(n)     = sum_{m!=n} V_i(nm) * V_j(mn) / (E_n - E_m)
%   |n^(1)>_i    = sum_{m!=n} |m> V_i(mn)/(E_n - E_m)
%
% IMPORTANT IMPLEMENTATION DETAILS:
%   • Denominators use (E_n - E_m) with indices arranged to match the formulae
%   • Projected matrices are explicitly Hermitianized to avoid numerical drift
%   • Degenerate PT: within each degenerate block, we rotate by diagonalizing V_i
%     Second-order sums exclude intra-block couplings; for cross-terms we return
%     effective block tensors if calcTensor=true
%
% EXAMPLES:
%   % Single perturbation
%   results = perturb2nd(H0, V);
%   
%   % Vector perturbation (e.g., magnetic field)
%   results = perturb2nd(H0, {Hx, Hy, Hz});
%   
%   % Along specific direction
%   results = perturb2nd(H0, {Hx, Hy, Hz}, 'direction', [1, 0, 1]);
%   
%   % With options
%   results = perturb2nd(H0, H_pert, 'threshold', 1e-10, 'verbose', false);
%
% SEE ALSO: eig

% Parse input arguments
p = inputParser;
p.FunctionName = 'perturb2nd';

addRequired(p, 'H0', @(x) isnumeric(x) && ismatrix(x) && size(x,1)==size(x,2));
addRequired(p, 'H_pert', @(x) (isnumeric(x) && ismatrix(x)) || iscell(x));

addParameter(p, 'threshold', 1e-12, @(x) isscalar(x) && x>=0);
addParameter(p, 'verbose', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'calcTensor', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'calcWaveFunc', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'direction', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));

parse(p, H0, H_pert, varargin{:});
opts = p.Results;

% ---- Normalize inputs ---------------------------------------------------
N = size(H0,1);

% Wrap single perturbation in a cell
if ~iscell(H_pert)
    Hlist = {H_pert};
else
    Hlist = H_pert(:).';
end
k = numel(Hlist);

% If a direction is provided, normalize its length and collapse to one component
useDirection = ~isempty(opts.direction);
if useDirection
    u = opts.direction(:);
    if numel(u) ~= k
        error('perturb2nd:badDirection', ...
            'Length of ''direction'' (%d) must match number of perturbations (%d).', numel(u), k);
    end
    if all(u==0)
        error('perturb2nd:badDirection', 'Direction vector must be non-zero.');
    end
    u = u / norm(u);
end

% ---- Hermitize inputs to be safe ---------------------------------------
H0 = hermitian_part(H0);
for i = 1:k
    Hlist{i} = hermitian_part(Hlist{i});
end

if opts.verbose
    fprintf('=== Second-Order Perturbation Theory ===\n');
    fprintf('System size: %d x %d\n', N, N);
    fprintf('Perturbation components: %d\n', k);
end

% ---- Solve unperturbed problem -----------------------------------------
if opts.verbose, fprintf('Computing unperturbed eigenvalues and eigenvectors...\n'); end
[W0, D0] = eig(H0);              % columns of W0 are eigenvectors
E0 = real(diag(D0));             % eigenvalues (may not be sorted)

% Sort energies ascending for stable grouping
[results.eigE0, perm] = sort(E0, 'ascend');
results.eigW0 = W0(:, perm);

% Project perturbations into eigenbasis of H0
results.Hproj = cell(1,k);
for i = 1:k
    H = results.eigW0' * Hlist{i} * results.eigW0;
    results.Hproj{i} = hermitian_part(H);  % enforce Hermiticity numerically
end

% ---- Detect degeneracies ------------------------------------------------
E = results.eigE0(:);
tol = opts.threshold;            % absolute energy tolerance
[blocks, hasDeg] = group_degenerate_levels(E, tol);

results.hasDegeneracy = hasDeg;
results.blocks = blocks;

if hasDeg && opts.verbose
    fprintf('  Note: Detected %d degenerate/near-degenerate block(s) (tol = %g).\n', numel(blocks), tol);
end

% ---- Combine along direction if requested ------------------------------
Hproj_eff = results.Hproj;
if useDirection
    % Collapse to a single effective perturbation in the eigenbasis
    H_eff = zeros(N,N);
    for i = 1:k
        H_eff = H_eff + u(i) * results.Hproj{i};
    end
    Hproj_eff = { (H_eff + H_eff')/2 };
    k_eff = 1;
    calcTensor = false;
else
    k_eff = k;
    calcTensor = opts.calcTensor;
end

% ---- Compute first/second-order corrections ----------------------------
if opts.verbose, fprintf('Computing perturbation corrections...\n'); end
[dE1, dE2, dPsi1, aux] = calcPertCorrections(E, results.eigW0, Hproj_eff, ...
                                             blocks, hasDeg, tol, opts.verbose, calcTensor, opts.calcWaveFunc);

% ---- Assemble outputs ---------------------------------------------------
results.dE1 = dE1;
results.dPsi1 = dPsi1;

if calcTensor
    % Full tensor K_ij per level (only if non-degenerate)
    results.dE2 = dE2;                    % [N x k x k]
    if hasDeg
        % Auxiliary block tensors are returned separately
        results.blockK = aux.blockK;
    else
        results.blockK = struct([]);
    end
else
    % Diagonal-only second order (per component or along direction)
    results.dE2 = dE2;                    % [N x size(dE1,2)]
    results.blockK = struct([]);
end

% Store parameters (echo)
results.parameters = opts;
if useDirection
    results.parameters.direction = u;
end

if opts.verbose
    fprintf('Perturbation theory calculation completed.\n');
    fprintf('=========================================\n');
end

end

%% Local Functions

function [blocks, hasDeg] = group_degenerate_levels(E, tol)
% Group sorted energies into degenerate blocks using a threshold on neighbours.
N = numel(E);
blocks = cell(1, 0);
if N == 0
    hasDeg = false;
    return;
end

idx = 1;
while idx <= N
    j = idx;
    while j < N && abs(E(j+1) - E(j)) < tol
        j = j + 1;
    end
    blocks{end+1} = idx:j;
    idx = j + 1;
end

hasDeg = any(cellfun(@(g) numel(g) > 1, blocks));
end
function [dE1, dE2, dPsi1, aux] = calcPertCorrections(E, W0, Hmats, ...
                                                      blocks, hasDeg, tol, verbose, calcTensor, calcWaveFunc)
% Compute first- and second-order corrections using projected matrices Hmats.
% E is [N x 1], W0 is [N x N], Hmats is {k_eff}[N x N]

N = numel(E);
k_eff = numel(Hmats);

% Hermitize once to avoid repeated numerical drift
Hherm = cellfun(@hermitian_part, Hmats, 'UniformOutput', false);

En_minus_Em = E(:) - E(:).';          % (n,m) entry = E_n - E_m
Em_minus_En = -En_minus_Em;           % (m,n) entry = E_m - E_n

offDiagMask = ~eye(N,'logical');

validNM = offDiagMask & (abs(En_minus_Em) > tol);
validMN = validNM.';

sameBlockNM = false(N);
if hasDeg
    for g = 1:numel(blocks)
        G = blocks{g};
        sameBlockNM(G,G) = true;
    end
end
sameBlockMN = sameBlockNM.';

maskE2 = validNM;
maskPsi = validMN;
if hasDeg
    maskE2 = maskE2 & ~sameBlockNM;
    maskPsi = maskPsi & ~sameBlockMN;
end

dE1 = zeros(N, k_eff);
if calcTensor
    dE2 = zeros(N, k_eff, k_eff);
else
    dE2 = zeros(N, k_eff);
end

if calcWaveFunc
    dPsi1 = cell(1, k_eff);
else
    dPsi1 = cell(1, 0);
end

aux = struct();
aux.blockK = struct([]);

for i = 1:k_eff
    Hi = Hherm{i};
    Ui = eye(N);
    diagHi = real(diag(Hi));
    if hasDeg
        for g = 1:numel(blocks)
            G = blocks{g};
            Hig = hermitian_part(Hi(G,G));
            [Ug, Dg] = eig(Hig);
            [dvals, ord] = sort(real(diag(Dg)), 'ascend');
            Ug = Ug(:, ord);
            Ui(G,G) = Ug;
            diagHi(G) = dvals;
        end
    end
    dE1(:,i) = diagHi;
    Hi_rot = Ui' * Hi * Ui;

    if ~calcTensor
        contrib = zeros(N);
        contrib(maskE2) = abs(Hi_rot(maskE2)).^2 ./ En_minus_Em(maskE2);
        dE2(:,i) = real(sum(contrib, 2));
    end

    if calcWaveFunc
        coeffs = zeros(N);
        coeffs(maskPsi) = Hi_rot(maskPsi) ./ Em_minus_En(maskPsi);
        Wrot = W0 * Ui;
        dPsi1{i} = Wrot * coeffs;
    end
end

if calcTensor
    if ~hasDeg
        for i = 1:k_eff
            Hi = Hherm{i};
            for j = 1:k_eff
                Hj = Hherm{j};
                HjT = Hj.'; % plain transpose to get (m,n)
                contrib = zeros(N);
                contrib(validNM) = Hi(validNM) .* HjT(validNM) ./ En_minus_Em(validNM);
                dE2(:,i,j) = real(sum(contrib, 2));
            end
        end
    else
        aux.blockK = repmat(struct('idx',[],'U',[],'K',[]), 1, numel(blocks));
        for g = 1:numel(blocks)
            G = blocks{g};
            notG = true(1,N); notG(G) = false;

            Hsum = zeros(numel(G));
            for i = 1:k_eff
                Hgg = Hherm{i}(G,G);
                Hsum = Hsum + Hgg*Hgg;
            end
            [Ublk, ~] = eig(hermitian_part(Hsum));

            HiG_notG = cell(1,k_eff);
            HjnotG_G = cell(1,k_eff);
            for i = 1:k_eff
                Hi = Hherm{i};
                HiG_notG{i} = Ublk' * Hi(G, notG);
                HjnotG_G{i} = Hi(notG, G) * Ublk;
            end

            Eblock = mean(E(G));
            denom = (Eblock - E(notG)).';

            Kcell = cell(k_eff, k_eff);
            for i = 1:k_eff
                weightedHi = HiG_notG{i} ./ denom;
                for j = 1:k_eff
                    Kblock = weightedHi * HjnotG_G{j};
                    if i == j
                        Kblock = hermitian_part(Kblock);
                    end
                    Kcell{i,j} = Kblock;
                end
            end

            aux.blockK(g).idx = G;
            aux.blockK(g).U   = Ublk;
            aux.blockK(g).K   = Kcell;
        end

        if verbose
            fprintf(['  calcTensor=true requested but degeneracies exist.\n', ...
                     '  Returned block-level tensors in results.blockK (basis in .U).\n']);
        end
    end
end
end

function Hh = hermitian_part(H)
% Return the Hermitian part of a numeric matrix.
Hh = (H + H')/2;
end





