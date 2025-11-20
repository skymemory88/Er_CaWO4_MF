function [en, ham_E, basis, jx, jy, jz] = Ising_proj(const, ion, params)
temp = params.temp;
Bfield = params.field;
gE = -ion.gLande(ion.idx) * const.muB * const.J2meV; % Zeeman interaction prefactor

en = double.empty(2, length(temp), size(Bfield,2), 0); % eigen-energy container
ham_E = double.empty(2, 2, length(temp), size(Bfield,2), 0); % eigen-energy container
basis = double.empty(2*ion.J(ion.idx)+1, 2, length(temp), size(Bfield,2), 0); % eigen-state container
jx = double.empty(1, length(temp), size(Bfield,2), 0); % <Jx> container
jy = double.empty(1, length(temp), size(Bfield,2), 0); % <Jy> container
jz = double.empty(1, length(temp), size(Bfield,2), 0); % <Jz> container
for ff = 1:size(Bfield,2)
    % construct the Hamiltonian
    Hz = gE*(Bfield(1,ff)*ion.Jx + Bfield(2,ff)*ion.Jy + Bfield(3,ff)*ion.Jz); % Electronic Zeeman interaction [meV]
    % Ham = ion.Hcf; % Crystal electrostatic field terms only
    % Ham = ion.Hcf + Hz; % Crystal field + Zeeman
    Ham = ion.Hcf + ion.h4 + Hz; % include h4 anisotropy term

    % Diagonalize
    [eigen, ee] = eig(Ham); % diagonalization by unitary transformation
    [ee, ~] = sort(real(diag(ee))); % Take only the real part of the eigen-energy to form a diaganol matrix
    ee = ee - min(ee);

    eigen = eigen(:,1:2); % truncate to Ising subspace

    % bias spin up for degenerate states (spontaneous symmetry breaking)
    spz = eigen' * ion.Jz * eigen; % truncate Jz operator
    [uni, spz] = eig(spz); % second rotation to diagonalize Jz
    if spz(1) < 0 % ensure consistent spontaneous symmetry breaking
        uni = flip(uni,2);
    end

    % fix the relative phase between the two basis vectors
    phi1 = angle(uni(1,1));
    uni(:,1) = uni(:,1) * exp(-1i*phi1);
    phi2 = angle(uni(1,2));
    uni(:,2) = uni(:,2) * exp(-1i*phi2);
    eigen = eigen * uni;

    for tt = 1:length(temp)
        beta = 1/const.kB/const.J2meV/temp(tt);
        if temp(tt) ~= 0
            Z = exp(-ee(1:2)* beta)/sum(exp(-ee(1:2) * beta)); % Truncate the partition function
            jx(:,tt,ff,1) = real( diag(eigen' * ion.Jx * eigen) )' * Z; % <Jx>
            jy(:,tt,ff,1) = real( diag(eigen' * ion.Jy * eigen) )' * Z; % <Jy>
            jz(:,tt,ff,1) = real( diag(eigen' * ion.Jz * eigen) )' * Z; % <Jz>
        else
            jx(:,tt,ff,1) = real( diag(eigen(:,1)' * ion.Jx * eigen(:,1)) )'; % <Jx>
            jy(:,tt,ff,1) = real( diag(eigen(:,1)' * ion.Jy * eigen(:,1)) )'; % <Jy>
            jz(:,tt,ff,1) = real( diag(eigen(:,1)' * ion.Jz * eigen(:,1)) )'; % <Jz>
        end
        ham_E(:,:,tt,ff,1) = eigen' * Ham * eigen; % matrix form in its eigen basis (diagonal)
        en(:,tt,ff,1) = ee(1:2); % ground state energy
        basis(:,:,tt,ff,1) = eigen; % truncate the basis to form an Ising basis
    end
end
% hamI = repmat(hamI,1,1,length(temp),1);
% en = repmat(en, 1, length(temp), 1);
% basis = repmat(basis,1,1,length(temp),1);
end