function [proj, Jproj] = projectDoublet(Hcf, Jx, Jy, Jz)
%PROJECTDOUBLET Project crystal-field operators onto lowest Kramers doublet
%   Returns the unitary that spans the lowest-energy doublet of Hcf and the
%   corresponding projected J operators expressed in that basis.

    % Enforce Hermiticity against numerical noise
    Hcf = (Hcf + Hcf') / 2;

    % Diagonalise and retain the lowest Kramers pair
    [eigVec, eigVal] = eig(Hcf);
    [evals, order] = sort(real(diag(eigVal)), 'ascend');
    eigVec = eigVec(:, order);
    evals = evals(:);
    subspace = eigVec(:, 1:2);

    % Diagonalise Jz within the doublet to obtain a reproducible gauge
    JzSub = (subspace' * Jz * subspace);
    JzSub = (JzSub + JzSub') / 2;
    [rot, Djz] = eig(JzSub);
    [~, rotOrder] = sort(real(diag(Djz)), 'descend');
    rot = rot(:, rotOrder);
    Djz = Djz(rotOrder, rotOrder);
    subspace = subspace * rot;

    % Fix arbitrary phases by making the first significant component real
    for ii = 1:2
        vec = subspace(:, ii);
        idx = find(abs(vec) > 1e-8, 1, 'first');
        if ~isempty(idx)
            phase = angle(vec(idx));
            subspace(:, ii) = vec * exp(-1i * phase);
        end
    end

    % Project J operators and re-symmetrise
    JprojX = subspace' * Jx * subspace;
    JprojY = subspace' * Jy * subspace;
    JprojZ = subspace' * Jz * subspace;
    JprojX = (JprojX + JprojX') / 2;
    JprojY = (JprojY + JprojY') / 2;
    JprojZ = (JprojZ + JprojZ') / 2;

    proj = struct();
    proj.states = subspace;
    proj.energies = evals(1:2);

    Jproj = struct();
    Jproj.Jx = JprojX;
    Jproj.Jy = JprojY;
    Jproj.Jz = JprojZ;
end
