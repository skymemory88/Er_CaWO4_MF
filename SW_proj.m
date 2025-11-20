function [en, ham_E, basis, jx, jy, jz] = SW_proj(const, ion, params)
% Schrieffer–Wolff (2nd order) truncation onto the lowest Kramers doublet
% Inputs/outputs and naming mirror your Ising_proj().

temp   = params.temp;
Bfield = params.field;                             % 3 x nFields
gE     = -ion.gLande(ion.idx) * const.muB * const.J2meV;  % [meV/T]

N   = 2*ion.J(ion.idx) + 1;                        % full Hilbert dimension
nT  = length(temp);
nB  = size(Bfield,2);

% Preallocate like your function (last dim kept as 1)
en    = zeros(2, nT, nB, 1);
ham_E = zeros(2, 2, nT, nB, 1);
basis = zeros(N, 2, nT, nB, 1);
jx    = zeros(1, nT, nB, 1);
jy    = zeros(1, nT, nB, 1);
jz    = zeros(1, nT, nB, 1);

% Convenient handles
Jx = ion.Jx; Jy = ion.Jy; Jz = ion.Jz;
H0 = ion.Hcf;                                    % unperturbed Hamiltonian

% If ion.h4 not provided, treat as zero
if isfield(ion,'h4') && ~isempty(ion.h4)
    H4 = ion.h4;
else
    H4 = zeros(N);
end

% Diagonalize H0 once (field-independent)
[U0, D0] = eig(H0); 
E0 = real(diag(D0)); 
[E0, ord] = sort(E0,'ascend'); 
U0 = U0(:, ord);

% Indices for P (lowest doublet) and Q (rest)
pIdx = 1:2;
qIdx = 3:N;

% Energies split
EP = E0(pIdx);                        % 2x1
EQ = E0(qIdx);                        % (N-2)x1

% Tolerance for denominators (avoid division by zero in pathological cases)
reg = 1e-12;

for ff = 1:nB
    % Zeeman term in the J basis
    Hz = gE * (Bfield(1,ff)*Jx + Bfield(2,ff)*Jy + Bfield(3,ff)*Jz);

    % Perturbation V in J basis (treat H4 + Hz as small vs H0)
    V = H4 + Hz;

    % Transform V to H0-eigenbasis
    V0 = U0' * V * U0;

    % Blocks of V0
    Vpp = V0(pIdx, pIdx);        % 2x2
    % Vpq = V0(pIdx, qIdx);        % 2x(N-2)
    Vqp = V0(qIdx, pIdx);        % (N-2)x2
    % Vqq = V0(qIdx, qIdx);      % not needed at 2nd-order symmetric SW

    % ----- 2nd-order degenerate SW effective Hamiltonian in P -----
    % Heff = E_P (diag) + V_pp + ΔH^(2)
    Heff = diag(EP) + Vpp;

    % ΔH^(2)_{ab} = - sum_r V_{a r} V_{r b} * 0.5*( 1/(E_a-E_r) + 1/(E_b-E_r) )
    DeltaH2 = zeros(2,2);
    for a = 1:2
        for b = 1:2
            denom_a = EP(a) - EQ;     % (N-2)x1
            denom_b = EP(b) - EQ;     % (N-2)x1
            % regularize tiny denominators (should not happen if CF gaps are finite)
            denom_a(abs(denom_a) < reg) = sign(real(denom_a(abs(denom_a) < reg))) * reg;
            denom_b(abs(denom_b) < reg) = sign(real(denom_b(abs(denom_b) < reg))) * reg;

            symDen = 0.5*(1./denom_a + 1./denom_b);  % (N-2)x1
            vec_ar = V0(pIdx(a), qIdx);              % 1x(N-2)
            vec_rb = V0(qIdx, pIdx(b));              % (N-2)x1
            % elementwise sum over r
            DeltaH2(a,b) = - sum( (vec_ar.' .* vec_rb) .* symDen );
        end
    end
    % enforce Hermiticity against tiny numerical drift
    DeltaH2 = (DeltaH2 + DeltaH2')/2;

    Heff = Heff + DeltaH2;

    % Diagonalize 2x2 effective Hamiltonian
    [Up, Deff] = eig(Heff);
    Eeff = real(diag(Deff));
    [Eeff, ord2] = sort(Eeff,'ascend');
    Up = Up(:, ord2);

    % ----- Build approximate full-space eigenvectors via SW mixing -----
    % S_QP(r,a) = V_{r a}/(E_a-E_r)   (in H0 basis)
    Den = (ones(length(EQ),1)*EP.') - (EQ*ones(1,2));     % (N-2)x2
    Den(abs(Den) < reg) = sign(real(Den(abs(Den) < reg))) * reg;
    Sqp = Vqp ./ Den;                                      % (N-2)x2

    % cP are coefficients of P-sector eigenvectors (columns of Up)
    cP = Up;                           % 2x2
    cQ = - Sqp * cP;                   % (N-2)x2

    % Full vector in H0 basis then back to J basis
    psi_H0 = [cP; cQ];                 % Nx2 (in H0 basis)
    eigen  = U0 * psi_H0;              % Nx2 (back in J basis)

    % Normalize columns
    for k = 1:2
        nk = sqrt(real(eigen(:,k)'*eigen(:,k)));
        if nk > 0, eigen(:,k) = eigen(:,k)/nk; end
    end

    % (Optional but matches your convention) Choose "spin-up/down" by Jz
    spz = eigen' * Jz * eigen;
    [uni, spz_vals] = eig(spz);
    % enforce consistent ordering (spz(1)<0 like your code)
    if real(spz_vals(1,1)) < 0
        uni = fliplr(uni);
    end
    % fix relative phases to make first component real and positive
    phi1 = angle(uni(1,1)); uni(:,1) = uni(:,1)*exp(-1i*phi1);
    phi2 = angle(uni(1,2)); uni(:,2) = uni(:,2)*exp(-1i*phi2);

    eigen = eigen * uni;

    % Energies & observables
    Eout = Eeff - min(Eeff);     % shift so ground is 0 meV, like your ee = ee - min(ee)

    for tt = 1:nT
        beta = 1/(const.kB*const.J2meV*temp(tt));
        if temp(tt) ~= 0
            w = exp(-Eout*beta);  w = w/sum(w);       % 2-level partition
            jx(:,tt,ff,1) = real(diag(eigen' * Jx * eigen)).' * w;
            jy(:,tt,ff,1) = real(diag(eigen' * Jy * eigen)).' * w;
            jz(:,tt,ff,1) = real(diag(eigen' * Jz * eigen)).' * w;
        else
            jx(:,tt,ff,1) = real(eigen(:,1)' * Jx * eigen(:,1));
            jy(:,tt,ff,1) = real(eigen(:,1)' * Jy * eigen(:,1));
            jz(:,tt,ff,1) = real(eigen(:,1)' * Jz * eigen(:,1));
        end
        ham_E(:,:,tt,ff,1) = diag(Eeff(ord2));   % diagonal in Heff eigenbasis
        en(:,tt,ff,1)      = Eout;
        basis(:,:,tt,ff,1) = eigen;              % Nx2 truncated basis in J basis
    end

    % (Optional sanity check): SW validity via Q-weight
    % mix_w = sum(abs(cQ).^2,1);  % should be << 1
    % fprintf('Field #%d: Q-weight ~ [%.3g, %.3g]\n', ff, mix_w(1), mix_w(2));

end
end
