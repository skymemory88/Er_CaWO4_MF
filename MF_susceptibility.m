% Magnetic susceptibility calculation
function results = MF_susceptibility(varargin)
% Handle both new and legacy calling conventions
if nargin == 7 && isstruct(varargin{1}) && isfield(varargin{1}, 'freq')
    % Distinguish by third argument: struct (params) vs scalar (J)
    if isstruct(varargin{3})
        % New signature: MF_susceptibility(Opts, const, params, eigE, eigW, Bf, spinOps)
        [Opts, const, params, eigE, eigW, Bf, spinOps] = varargin{:};
    else
        % Legacy signature: MF_susceptibility(Options, const, J, I, eigenE, eigenW, Bfield)
        [Options, const, J, I, eigenE, eigenW, Bfield] = varargin{:};
        % Convert to new format
        Opts = Options;
        if ~isfield(Opts, 'parallel')
            Opts.parallel = false; % Default to serial for backward compatibility
        end
        params.J = J;
        params.I = I;
        eigE = eigenE;
        eigW = eigenW;
        Bf = Bfield;
        % Generate spin operators on the fly
        [~,~,~,~,~,~,Jxh,Jyh,Jzh,Ixh,Iyh,Izh] = spin_operators(J, I);
        spinOps = struct('Jxh', Jxh, 'Jyh', Jyh, 'Jzh', Jzh, 'Ixh', Ixh, 'Iyh', Iyh, 'Izh', Izh);
    end
else
    error('Invalid number of arguments. Expected 7 arguments with Options/Opts as first argument.');
end
% MF_SUSCEPTIBILITY Calculate magnetic susceptibility tensor
%
% NEW SIGNATURE:
%   results = MF_susceptibility(Opts, const, params, eigE, eigW, Bf, spinOps)
%
% LEGACY SIGNATURE (for backward compatibility):
%   results = MF_susceptibility(Options, const, J, I, eigenE, eigenW, Bfield)
%
% Inputs (new signature):
%   Opts    - Options structure with fields:
%             .freq     - Frequency array [GHz]  
%             .temp     - Temperature (scalar or array)
%             .gamma    - Broadening parameter [meV]
%             .parallel - Enable parallel processing (optional, default: false)
%   const   - Constants structure with .Gh2mV and .kB_meV
%   params  - Parameters structure with .J and .I
%   eigE    - Eigenvalues (ns x nb)
%   eigW    - Eigenvectors (ns x ns x nb)
%   Bf      - Magnetic field array (3 x nb)
%   spinOps - Pre-computed spin operators structure
%
% Inputs (legacy signature):
%   Options - Options structure (same as Opts above, .parallel optional)
%   const   - Constants structure
%   J       - Total angular momentum quantum number
%   I       - Nuclear spin quantum number
%   eigenE  - Eigenvalues (ns x nb) 
%   eigenW  - Eigenvectors (ns x ns x nb)
%   Bfield  - Magnetic field array (3 x nb)
%
% Output:
%   results - Structure with susceptibility tensor and metadata

fprintf('\nCalculating magnetic susceptibility (optimized)...\n');

% Setup grids
freq = Opts.freq(:)'; % [GHz]
omega = freq * const.Gh2mV; % [meV]

% Get dimensions
[~, nb] = size(eigE);
nf = length(omega);

% Setup temperatures
if isscalar(Opts.temp)
    temps = repmat(Opts.temp, 1, nb);
else
    temps = Opts.temp(:)';
    assert(length(temps) == nb, 'Temperature array size mismatch');
end

% Use magnetic moment operators (electronic + nuclear)
if params.I > 0
    gN_eff = 0.1;  % effective nuclear g-factor weighting
    JhT_x = spinOps.Jxh + gN_eff * spinOps.Ixh;
    JhT_y = spinOps.Jyh + gN_eff * spinOps.Iyh;
    JhT_z = spinOps.Jzh + gN_eff * spinOps.Izh;
else
    JhT_x = spinOps.Jxh;
    JhT_y = spinOps.Jyh;
    JhT_z = spinOps.Jzh;
end
kB_meV = const.kB_meV;
gamma = Opts.gamma; % [meV]

% Initialize output array
chi = complex(zeros(3, 3, nf, nb));
if Opts.parallel
    fprintf('Using parallel processing for %d field points...\n', nb);
    parfor ii = 1:nb
        v_loc = eigW(:, :, ii);
        en_loc = eigE(:, ii);
        temp_loc = temps(ii);

        chi_loc = chi_tensor(v_loc, en_loc, temp_loc, omega, ...
            JhT_x, JhT_y, JhT_z, kB_meV, gamma, nf);

        chi(:, :, :, ii) = chi_loc;
    end
else
    fprintf('Using optimized serial processing...\n');
    progInt = max(1, floor(nb / 10));
    for ii = 1:nb
        if mod(ii, progInt) == 0
            fprintf('Processing field point %d/%d (%.1f%%)\n', ii, nb, 100*ii/nb);
        end

        v = eigW(:, :, ii);
        en = eigE(:, ii);
        temp = temps(ii);

        chi(:, :, :, ii) = chi_tensor(v, en, temp, omega, ...
            JhT_x, JhT_y, JhT_z, kB_meV, gamma, nf);
    end
end

results.chi = chi;
results.freq = freq;  % Store in GHz
results.omega = omega;  % Store in meV
results.temp = temps;
results.Bfield = Bf;  % Store full field

plotSusceptibility(chi, Bf, omega);
fprintf('Calculation completed.\n');
end

% Core susceptibility calculation - vectorized over frequencies
function chi_ij = chi_tensor(v, en, temp, omega, JhT_x, JhT_y, JhT_z, kB_meV, gamma, nf)
ns = length(en);

% Thermal population
if temp > 0
    beta = 1 / (kB_meV * temp);
    Z = sum(exp(-beta * en));
    zn = exp(-beta * en) / Z;
else % T=0 limit: only ground state populated
    zn = zeros(size(en));
    zn(1) = 1;
end

% Population difference matrix
[n, np] = meshgrid(zn, zn);
NN = n - np;

% Matrix elements in eigenstate basis
ttx = v' * JhT_x * v;
tty = v' * JhT_y * v;
ttz = v' * JhT_z * v;

% Energy differences
[ee, eep] = meshgrid(en, en);
eDiff = eep - ee;

% Vectorize over frequencies
eDiff3d = repmat(eDiff, [1, 1, nf]);
omega3d = reshape(omega, [1, 1, nf]);
omega3d = repmat(omega3d, [ns, ns, 1]);

% Denominator with broadening
EE_full = eDiff3d - omega3d;
deno = 1 ./ (EE_full - 1i * gamma);

% Expand matrices to 3D for vectorized calculation
NN3d = repmat(NN, [1, 1, nf]);
ttx3d = repmat(ttx, [1, 1, nf]);
tty3d = repmat(tty, [1, 1, nf]);
ttz3d = repmat(ttz, [1, 1, nf]);

% Calculate all susceptibility tensor components
chi_ij = complex(zeros(3, 3, nf));

% Diagonal components
chi_ij(1,1,:) = squeeze(sum(sum(ttx3d .* conj(ttx3d) .* NN3d .* deno, 1), 2));
chi_ij(2,2,:) = squeeze(sum(sum(tty3d .* conj(tty3d) .* NN3d .* deno, 1), 2));
chi_ij(3,3,:) = squeeze(sum(sum(ttz3d .* conj(ttz3d) .* NN3d .* deno, 1), 2));

% Off-diagonal components
chi_ij(1,2,:) = squeeze(sum(sum(ttx3d .* conj(tty3d) .* NN3d .* deno, 1), 2));
chi_ij(1,3,:) = squeeze(sum(sum(ttx3d .* conj(ttz3d) .* NN3d .* deno, 1), 2));
chi_ij(2,1,:) = squeeze(sum(sum(tty3d .* conj(ttx3d) .* NN3d .* deno, 1), 2));
chi_ij(2,3,:) = squeeze(sum(sum(tty3d .* conj(ttz3d) .* NN3d .* deno, 1), 2));
chi_ij(3,1,:) = squeeze(sum(sum(ttz3d .* conj(ttx3d) .* NN3d .* deno, 1), 2));
chi_ij(3,2,:) = squeeze(sum(sum(ttz3d .* conj(tty3d) .* NN3d .* deno, 1), 2));
end

% Plotting the susceptibility
function plotSusceptibility(chi, Bf, omega)
% Calculate field magnitude with sign
B0 = vecnorm(Bf, 2, 1) .* sign(sum(Bf, 1));

% Create mesh for plotting
[Bmesh, Fmesh] = meshgrid(B0 * 1000, omega);  % Convert to mT

% Extract components
chiXXreal = real(squeeze(chi(1,1,:,:)));
chiXXimag = imag(squeeze(chi(1,1,:,:)));
chiZZreal = real(squeeze(chi(3,3,:,:)));
chiZZimag = imag(squeeze(chi(3,3,:,:)));

% Create figure
figure('Name', 'Magnetic Susceptibility');

% Create subplots
subplot(2,2,1);
formatSubplot(Bmesh, Fmesh, chiXXreal, 'Re[\chi_{xx}]');

subplot(2,2,2);
formatSubplot(Bmesh, Fmesh, chiXXimag, 'Im[\chi_{xx}]');

subplot(2,2,3);
formatSubplot(Bmesh, Fmesh, chiZZreal, 'Re[\chi_{zz}]');

subplot(2,2,4);
formatSubplot(Bmesh, Fmesh, chiZZimag, 'Im[\chi_{zz}]');

colormap('jet');
sgtitle('Magnetic Susceptibility Tensor Components', 'FontSize', 14);
end

% Helper function for subplot formatting
function formatSubplot(Bmesh, Fmesh, data, titleStr)
hp = pcolor(Bmesh, Fmesh, data);
set(hp, 'EdgeColor', 'none');
title(titleStr, 'FontSize', 11);
xlabel('B (mT)');
ylabel('\omega (meV)');
colorbar;
axis tight;
end