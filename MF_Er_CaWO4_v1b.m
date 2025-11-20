%% Main Script
clearvars

% Setup options
Options.CEF = true; % use crystal field parameters
Options.nEn = 16; % number of eigenstates to include
Options.ndE = [1:15]; % specify which transition orders to calculate and plot
% Options.ndE = [4, 6, 8, 10, 11, 12, 13, 14, 15]; % specify which transition orders to calculate and plot
Options.dT_dB = false; % First derivatives of transitions
Options.zefoz = [-28E-6, 28E-6, 3.8E-6]; % ZEFOZ criteria (field range, criteria)
Options.d2T_dB2 = false; % Second derivatives of transitions
Options.chi = false; % calculate magnetic susceptibility
Options.freq = linspace(0.3, 0.35, 101); % frequency range for susceptibility [GHz]
Options.temp = 0.1; % temperature for susceptibility calculation [K]
Options.gamma = 1e-6; % spin linewidth [meV]

Options.line_styles = {'-', '--', '-.', ':'}; % line styles for different transitions
Options.colors = {'b', 'r', 'g', 'm', 'c', 'k', 'y'}; % colors for different transitions

B0 = linspace(-0.05, 0.05, 101); % [T]
% B0 = [0 0.2]; % [T]
theta = 0; % [rad] deviation angle from c-axis
phi = 0; % [rad] rotation angle within ab-plane, relative to a-axis
Bfield = [B0*sin(theta)*cos(phi); B0*sin(theta)*sin(phi); B0*cos(theta)];

% Initialize system parameters and operators
[const, params, ion] = initialization(Bfield, Options);
% En_max = (2*params.J+1)*(2*params.I+1); % maximal number of states

% Calculate eigenvalues and eigenvectors
[eigenE, eigenW, dEn] = eigenEnergy(Bfield, Options, const, params, ion);

% Calculate transitions and derivatives - Modified to return mode_info
[deltaEn, mode_info, zefoz_modes] = deltaE(Options, ion.J, ion.I, dEn, Bfield);

% Display summary of filtered transitions
if ~isempty(zefoz_modes) && isfield(zefoz_modes, 'mode_order')
    fprintf('\n=== ZEFOZ Filtered Transitions Summary ===\n');

    % Get field indices for the window (needed for statistics)
    B0 = vecnorm(Bfield,2,1) .* sign(sum(Bfield,1));
    B_range = sort(Options.zefoz(1:2));
    field_mask = B0 >= B_range(1) & B0 <= B_range(2);
    field_indices = find(field_mask);

    for i = 1:length(zefoz_modes.mode_order)
        fprintf('Transition %d: |%d> → |%d> (Order %d)\n', ...
            i, zefoz_modes.iState(i), ...
            zefoz_modes.fState(i), ...
            zefoz_modes.mode_order(i));

        % Check for curvature at window center
        window_center_idx = field_indices(round(length(field_indices)/2));
        if window_center_idx > 0 && window_center_idx <= length(zefoz_modes.d2T_dB2{i})
            d2_at_center = zefoz_modes.d2T_dB2{i}(window_center_idx);
            if ~isnan(d2_at_center)
                if abs(d2_at_center) < 1e-6
                    fprintf('  Curvature at window center: ~flat (d²E/dB² ≈ 0)\n');
                elseif d2_at_center > 0
                    fprintf('  Curvature at window center: concave up (d²E/dB² = %.3e)\n', d2_at_center);
                else
                    fprintf('  Curvature at window center: concave down (d²E/dB² = %.3e)\n', d2_at_center);
                end
            end
        end
    end
    fprintf('==========================================\n');
end

% Calculate susceptibility if requested
if Options.chi
    results = MF_chi(Options, const, ion.J, ion.I, eigenE, eigenW, Bfield);
end

function [const, params, ion] = initialization(Bfield, Options)
% define constants
const.hbar = 1.05457E-34; % Reduced Planck constant [J.s]
const.muB = 9.274e-24; % [J/T]
const.muN = 5.05078e-27; % [J/T]
const.kB = 1.3806e-23; % [J/K]
const.mu0 = 4e-7 * pi; % [H/m]
const.kB_meV = 8.61733e-2; % Boltzmann constant [meV/K]
const.J2meV = 6.24151e+21; % [mev/J]
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV; % [meV/GHz]

ion.I = 7/2; % Er nuclear spin;

% nuclear g-tensor
params.gN = [0.1618 0.1618 0.1618]; % nuclear g-factor
% params.gN = [0 0 0]; % nuclear g-factor
% params.gN = [-1.36503716E-1 -1.36503716E-1 -1.36503716E-1]; % nuclear g-factor

% quadrupolar interaction tensor
params.Q = [1.67 1.67 -3.34] / 1000 * const.Gh2mV; % [meV] nuclear quadrupler interaction strength
% params.Q = [-15.73 -15.73 -20.73] / 1000 * const.Gh2mV; % [meV] nuclear quadrupler interaction strength
% params.Q = [0 0 0]; % nuclear quadrupler interaction strength
% params.Q = [-2.55041373 -2.55041373 -5.28392734] / 1000 * const.Gh2mV; % nuclear quadrupler interaction strength

if Options.CEF == true % explicit CEF hamiltonian
    % Option 1: use CEF
    L = 6; S = 3/2; % Er
    ion.gLande = gLande(L,S);
    ion.J = L + S;
    ion.h4 = 0; % h4 anisotropy term
    % ion.h4 = 5.26e-3; % h4 anisotropy term (Neda's thesis)
    [~,~,~,~,~,~,ion.Jx,ion.Jy,ion.Jz,~,~,~] = spin_operators(ion.J, 0);
    
    % crystal electric field (CEF) parameters
    % B = [133.0,  -164.0,  186.0,   0,   -4.8,   281.8,  0];  % [cm^-1] J. Chem. Phys. 55(5), 1971
    % B = [168.2829   16.3804   86.5702    0.0109    0.1101   -0.9408    0.0010];  % [cm^-1] arXiv:2412.03948 (2025)-1
    % B = [137.2615    0.4310   -0.2058    0.0169    0.0011   -0.0048   -0.0005];  % [cm^-1] arXiv:2412.03948 (2025)-2
    % B = [ 145.3067    4.1297   22.5464    0.0153   -0.0198    0.4539   -0.0002];  % [cm^-1] arXiv:2412.03948 (2025)-3
    % B = [ 129.9484   -3.7035  -22.5857    0.0186    0.0118   -0.4555   -0.0011];  % [cm^-1] arXiv:2412.03948 (2025)-4
    
    B = [-979.402 1898.05 346.52 -272.893 -7.86205 -187.516 58.9021]; % 2025.11.02
    % B = [-330.74 2212.59 196.501 5.57413 -8.07197 -217.553 81.2414]; % 2025.11.10
    % B = [-753.279 796.633 -376.577 -133.463 -3.30013 -84.5397 -9.7149]; % (2025.10.13)
    % B = [119.6, -146.1, 187.6, 0, -6.5, 284.3, 0];  % in cm^{-1}
    % B = [567.02 -945.07 1008.10 0 -22.16  936.29  0.854]; % [cm-1] J. Chem. Phys. 53 (9), 1970
    % B = [220 -79.5 -593 -582 -2.545 -314 -294]; % [cm-1] 
    % B = [548.81 -942.41 1004.66 0 -17.12  946.84  1.343]; % [cm-1] J. Chem. Phys. 53 (9), 1970
    % B = [238 -90 -852 0 -0.6 -396 -75]; % [cm-1] J. Chem. Phys. 55, 5 (1971)
    % B = [640 -809 860 0 12  444  161]; % [cm-1] J. Chem. Phys. 54 314, (1971)
    % B = [732 -851 963 0 -9  555  105]; % [cm-1] J. Chem. Phys. 55, 2538 (1971)
    B = B .* 0.123983; % [cm-1] --> [meV]

    ion.Hcf = cf(ion.J, B, 0); % Build crystal field hamiltonian
    ion.idx = 1; % Er (legacy index)
    params.temp = Options.temp;
    ham_E = double.empty(2,2,size(Bfield,2),0);
    basis = double.empty(2*ion.J+1,2,size(Bfield,2),0);
    for ii = 1:size(Bfield,2)
        params.field = Bfield(:,ii);
        [~, ham_E(:,:,ii,1), basis(:,:,ii,1), ~, ~, ~] = SW_proj(const, ion, params);
        % [~, ham_E(:,:,ii,1), basis(:,:,ii,1), ~, ~, ~] = Ising_proj(const, ion, params);
    end
    % Map effective spin-1/2 hyperfine to the isotropic CEF scale:
    % A_eff = A0 * g_eff / gJ  =>  A0 = A_eff * gJ / g_eff
    % Target sets: A_eff = [-871.1, -871.1, -130.3] MHz, g_eff = [8.3, 8.3, 1.26],
    % gJ = gLande(6, 3/2) ≈ 1.20003 gives A0 ≈ -125.9 MHz
    A0_MHz = -125.9; % Scalar hyperfine constant consistent with target A_eff
    A = (A0_MHz/1000) * const.Gh2mV * [1 1 1];
    params.gE = ion.gLande * [1 1 1];
    params.basis = basis;
    ion.ham_E = ham_E;
else % effective spin-1/2 model
    ion.J = 1/2; % effective Kramer's spin
    
    % % Parameter set 1
    % A = -[873 873 130] / 1000 * const.Gh2mV; % [meV] hyperfine interaction
    % gE = [8.38 8.38 1.247]; % electronic g-factor
    
    % Parameter set 2
    A = -[871.1 871.1 130.3] / 1000 * const.Gh2mV; % [meV] hyperfine interaction
    gE = [8.3 8.3 1.26]; % electronic g-factor
    params.gE = gE;

    % Parameter set 3
    % A = [-871.436227 -871.436227 110.543961] / 1000 * const.Gh2mV; % [meV]
    % gE = [8.38 8.38 1.26234449]; % electronic g-factor
end
params.A = A;
end

function [eigenE, eigenW, dEn] = eigenEnergy(Bfield, Options, const, params, ion)
J = ion.J;
I = ion.I;
A = params.A;
Q = params.Q;
gE = params.gE;
gN = params.gN;

if Options.CEF == true
    % Nuclear spin operators
    Iz = diag(I:-1:-I);
    Ip = diag(sqrt((I - (I-1:-1:-I)) .* (I+1 + (I-1:-1:-I))), 1);
    Im = Ip';
    Ix = (Ip + Im)/2;
    Iy = (Ip - Im)/2i;

    Ixh = kron(eye(2),Ix);
    Iyh = kron(eye(2),Iy);
    Izh = kron(eye(2),Iz);
else
    En_max = (2*J+1)*(2*I+1); % maximal number of states
    if Options.nEn > En_max; Options.nEn = En_max; end
    [~,~,~,~,~,~,Jxh,Jyh,Jzh,Ixh,Iyh,Izh] = spin_operators(J, I);
end

% Initialize arrays
eigenE = zeros(Options.nEn, size(Bfield,2));
eigenW = zeros(Options.nEn, Options.nEn, size(Bfield,2));
dEn = zeros(Options.nEn, Options.nEn, size(Bfield,2));
for ii = 1:size(Bfield,2)
    if Options.CEF
        wav0 = params.basis(:,:,ii);
        jx = real(wav0' * ion.Jx * wav0);
        jy = real(wav0' * ion.Jy * wav0);
        jz = real(wav0' * ion.Jz * wav0);

        Jxh = kron(jx,eye(2*I+1));
        Jyh = kron(jy,eye(2*I+1));
        Jzh = kron(jz,eye(2*I+1));
        ham_E = squeeze(ion.ham_E(:,:,ii));
        ham_E = kron(ham_E,eye(2*I+1));
    else
        % electronic Zeeman interaction
        ham_E = -const.J2meV * const.muB *...
            (gE(1) * Jxh * Bfield(1,ii) + gE(2) * Jyh * Bfield(2,ii) + gE(3) * Jzh * Bfield(3,ii));
    end
    % nuclear Zeeman interaction
    ham_N = -const.J2meV * const.muN *...
        (gN(1) * Ixh * Bfield(1,ii) + gN(2) * Iyh * Bfield(2,ii) + gN(3) * Izh * Bfield(3,ii));
    % hyperfine interaction
    ham_hyp = A(1)*Ixh*Jxh + A(2)*Iyh*Jyh + A(3)*Izh*Jzh;
    % quadrupolar nuclear interaction 
    ham_QDP = Q(1)*Ixh*Ixh + Q(2)*Iyh*Iyh + Q(3)*Izh*Izh; % with quadrupolar interaction
 
    % total hamiltonian
    ham = ham_E + ham_N + ham_hyp + ham_QDP;
    % ham =  ham_E + ham_N + ham_hyp;

    [wv, ee] = eig(ham); 
    ee = real(diag(ee)); % [meV] Take only the real part of the eigen-energy to form a diagonal matrix
    [ee, n] = sort(ee); % sort the energy from lowest to the highest
    wv = wv(:,n);
    % ee = ee - min(ee); % [meV] Normalize the energy amplitude to the lowest eigen-energy
    eigenE(:,ii) = ee(1:Options.nEn); % [meV] 
    eigenW(:,:,ii) = wv(1:Options.nEn,1:Options.nEn); % sort the eigen-vectors in its basis accordingly

    % calculate all the possible inter-state transitions
    en = ee(1:Options.nEn) / const.Gh2mV; % [meV] --> [GHz]
    [Ex, Ey] = meshgrid(en,en);
    dEn(:,:,ii) = Ex - Ey; % [GHz] all possible transitions among eigenstates
end
end

function [deltaEn, mode_info, zefoz_modes] = deltaE(Options, J, I, dEn, Bfield)
En_max = (2*J+1)*(2*I+1);
if Options.nEn > En_max; Options.nEn = En_max; end

% Calculate transitions for specified orders
deltaEn = cell(length(Options.ndE), 1);
mode_info = cell(length(Options.ndE), 1);
for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < Options.nEn
        num_transitions = Options.nEn - trans_order;
        deltaEn{t_idx} = zeros(num_transitions, size(dEn,3));

        mode_info{t_idx}.iState = zeros(num_transitions, 1);
        mode_info{t_idx}.fStates = zeros(num_transitions, 1);

        for ii = 1:size(dEn,3)
            deltaEn{t_idx}(:, ii) = diag(squeeze(dEn(:,:,ii)), trans_order);
        end

        for trans_idx = 1:num_transitions
            mode_info{t_idx}.iState(trans_idx) = trans_idx - 1;
            mode_info{t_idx}.fStates(trans_idx) = trans_idx - 1 + trans_order;
        end
    end
end

fprintf('Calculated transitions for orders: %s\n', mat2str(Options.ndE));
fprintf('Total number of eigenstates: %d\n', Options.nEn);

% % Print transition information
% fprintf('\n=== Transition State Pairs ===\n');
% for t_idx = 1:length(Options.ndE)
%     trans_order = Options.ndE(t_idx);
%     if trans_order < Options.nEn && ~isempty(mode_info{t_idx})
%         fprintf('\nOrder %d transitions (|n> → |n+%d>):\n', trans_order, trans_order);
%         for trans_idx = 1:length(mode_info{t_idx}.iState)
%             fprintf('  Transition %d: |%d> → |%d>\n', trans_idx, ...
%                 mode_info{t_idx}.iState(trans_idx), ...
%                 mode_info{t_idx}.fStates(trans_idx));
%         end
%     end
% end
fprintf('==============================\n\n');

% Plot basic transitions with labels
plot_basic_transitions(Options, deltaEn, Bfield, En_max, mode_info);

% Calculate derivatives and ZEFOZ filtering if requested
zefoz_modes = [];
if Options.dT_dB
    % This will now include ZEFOZ filtering if Options.zefoz is defined
    zefoz_modes = dT_dB(Options, deltaEn, Bfield, En_max, mode_info);
end
end

function plot_basic_transitions(Options, deltaEn, Bfield, En_max, mode_info)
% Plotting with state labels
figure;
box on
hold on

plot_handles = [];
legend_labels = {};

% Extend line styles and colors if needed
num_transitions = length(Options.ndE);
if num_transitions > length(Options.line_styles)
    Options.line_styles = repmat(Options.line_styles, 1, ceil(num_transitions/length(Options.line_styles)));
end
if num_transitions > length(Options.colors)
    Options.colors = repmat(Options.colors, 1, ceil(num_transitions/length(Options.colors)));
end

B_plot = 1000 * vecnorm(Bfield,2,1) .* sign(sum(Bfield,1)); % [mT]
for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(deltaEn{t_idx})
        % Get the data for this transition order
        En_x = deltaEn{t_idx};

        if ~isempty(En_x)
            num_trans = size(En_x, 1);

            % Plot each transition with a unique color/style combination
            for trans_idx = 1:num_trans
                % Use different shades of the base color for transitions within same order
                color_shade = Options.colors{mod(t_idx-1, length(Options.colors))+1};
                if ischar(color_shade)
                    % Convert character color to RGB and adjust brightness
                    switch color_shade
                        case 'b', base_color = [0 0 1];
                        case 'r', base_color = [1 0 0];
                        case 'g', base_color = [0 1 0];
                        case 'm', base_color = [1 0 1];
                        case 'c', base_color = [0 1 1];
                        case 'k', base_color = [0 0 0];
                        case 'y', base_color = [1 1 0];
                        otherwise, base_color = [0 0 1];
                    end
                    % Adjust brightness for different transitions
                    brightness_factor = 0.7 + 0.3 * (trans_idx-1)/(max(num_trans-1,1));
                    plot_color = base_color * brightness_factor;
                    plot_color(plot_color > 1) = 1; % Cap at 1
                else
                    plot_color = color_shade;
                end

                % Determine line style - cycle through available styles for different transitions
                linStyl = Options.line_styles{mod(trans_idx-1, length(Options.line_styles))+1};

                p = plot(B_plot, En_x(trans_idx, :), linStyl, ...
                    'Color', plot_color, 'LineWidth', 1.5);

                % Create label with specific state information
                iState = mode_info{t_idx}.iState(trans_idx);
                fState = mode_info{t_idx}.fStates(trans_idx);
                label = sprintf('|%d⟩ → |%d⟩', iState, fState);

                plot_handles = [plot_handles, p];
                legend_labels{end+1} = label;
            end
        end
    end
end

% Customize plot
xlabel('Magnetic Field (mT) $\|$ c','interpreter','latex')
ylabel('Energy gap (GHz)')
set(gca,'fontsize',14)

% Add legend with specific state labels
if ~isempty(plot_handles)
    legend(plot_handles, legend_labels, 'Location', 'bestoutside', 'FontSize', 10)
end

grid on
title('Eigenstate Transitions vs Magnetic Field (with State Labels)')
end

function zefoz_modes = dT_dB(Options, deltaEn, Bfield, En_max, mode_info)
% Calculate first and second derivatives of transitions
fprintf('\nCalculating transition derivatives...\n');

B0 = vecnorm(Bfield,2,1) .* sign(sum(Bfield,1)); % [T]
dT = cell(length(Options.ndE), 1);
d2T = cell(length(Options.ndE), 1);

% Check if grid is uniform (it should be since linspace is used)
dB_grid = diff(B0);
is_uniform = std(dB_grid) < 1e-10 * mean(dB_grid);
if ~is_uniform
    warning('Non-uniform magnetic field grid detected. Derivative accuracy may be reduced.');
end

for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(deltaEn{t_idx})
        data_to_process = deltaEn{t_idx};
        num_transitions = size(data_to_process, 1);
        num_fields = size(data_to_process, 2);

        % Initialize derivative arrays
        dT{t_idx} = zeros(num_transitions, num_fields);
        if Options.d2T_dB2
            d2T{t_idx} = zeros(num_transitions, num_fields);
        end

        % Calculate derivatives for each transition
        for trans_idx = 1:num_transitions
            transition_freq = data_to_process(trans_idx, :);

            % First derivative calculation
            derivative = zeros(1, num_fields);

            if num_fields > 2
                % For uniform grid, use standard finite difference formulas
                if is_uniform
                    dB = B0(2) - B0(1); % [T]

                    % Forward difference for first point (2nd order accurate)
                    derivative(1) = (-3*transition_freq(1) + 4*transition_freq(2) - transition_freq(3)) / (2 * dB * 1000);

                    % Central difference for interior points (2nd order accurate)
                    for ii = 2:num_fields-1
                        derivative(ii) = (transition_freq(ii+1) - transition_freq(ii-1)) / (2 * dB * 1000);
                    end

                    % Backward difference for last point (2nd order accurate)
                    derivative(end) = (transition_freq(end-2) - 4*transition_freq(end-1) + 3*transition_freq(end)) / (2 * dB * 1000);

                else
                    % For non-uniform grid, use point-by-point calculation
                    % Forward difference for first point
                    dB = B0(2) - B0(1);
                    derivative(1) = (transition_freq(2) - transition_freq(1)) / (dB * 1000);

                    % Central difference for middle points
                    for ii = 2:num_fields-1
                        dB_forward = B0(ii+1) - B0(ii);
                        dB_backward = B0(ii) - B0(ii-1);

                        % Use weighted central difference for non-uniform grid
                        weight_forward = dB_backward / (dB_forward + dB_backward);
                        weight_backward = dB_forward / (dB_forward + dB_backward);

                        derivative(ii) = (weight_forward * (transition_freq(ii+1) - transition_freq(ii)) / dB_forward + ...
                            weight_backward * (transition_freq(ii) - transition_freq(ii-1)) / dB_backward) / 1000;
                    end

                    % Backward difference for last point
                    dB = B0(end) - B0(end-1);
                    derivative(end) = (transition_freq(end) - transition_freq(end-1)) / (dB * 1000);
                end

            elseif num_fields == 2
                % Only two points - use simple difference
                dB = B0(2) - B0(1);
                derivative = (transition_freq(2) - transition_freq(1)) / (dB * 1000) * ones(1, 2);
            else
                % Single point - derivative is undefined
                derivative = NaN;
            end

            dT{t_idx}(trans_idx, :) = derivative;

            % Second derivative calculation
            if Options.d2T_dB2
                second_derivative = zeros(1, num_fields);

                if num_fields > 3
                    if is_uniform
                        dB = B0(2) - B0(1); % [T]
                        dB2 = (dB * 1000)^2; % Convert to (mT)^2

                        % Forward difference for first point (2nd order accurate)
                        second_derivative(1) = (2*transition_freq(1) - 5*transition_freq(2) + 4*transition_freq(3) - transition_freq(4)) / dB2;

                        % Central difference for interior points (2nd order accurate)
                        for ii = 2:num_fields-1
                            second_derivative(ii) = (transition_freq(ii+1) - 2*transition_freq(ii) + transition_freq(ii-1)) / dB2;
                        end

                        % Backward difference for last point (2nd order accurate)
                        second_derivative(end) = (-transition_freq(end-3) + 4*transition_freq(end-2) - 5*transition_freq(end-1) + 2*transition_freq(end)) / dB2;

                    else
                        % For non-uniform grid, use simpler approximations
                        % Forward difference for first two points
                        dB1 = (B0(2) - B0(1)) * 1000; % [mT]
                        dB2 = (B0(3) - B0(2)) * 1000; % [mT]
                        f0 = transition_freq(1);
                        f1 = transition_freq(2);
                        f2 = transition_freq(3);

                        % Approximate second derivative at first point
                        second_derivative(1) = 2 * ((f2 - f1)/dB2 - (f1 - f0)/dB1) / (dB1 + dB2);

                        % Central difference for middle points
                        for ii = 2:num_fields-1
                            dB_back = (B0(ii) - B0(ii-1)) * 1000; % [mT]
                            dB_forw = (B0(ii+1) - B0(ii)) * 1000; % [mT]

                            deriv_back = (transition_freq(ii) - transition_freq(ii-1)) / dB_back;
                            deriv_forw = (transition_freq(ii+1) - transition_freq(ii)) / dB_forw;

                            second_derivative(ii) = 2 * (deriv_forw - deriv_back) / (dB_back + dB_forw);
                        end

                        % Backward difference for last point
                        dB1 = (B0(end-1) - B0(end-2)) * 1000; % [mT]
                        dB2 = (B0(end) - B0(end-1)) * 1000; % [mT]
                        f0 = transition_freq(end-2);
                        f1 = transition_freq(end-1);
                        f2 = transition_freq(end);

                        second_derivative(end) = 2 * ((f2 - f1)/dB2 - (f1 - f0)/dB1) / (dB1 + dB2);
                    end

                elseif num_fields == 3
                    % Exactly 3 points - can calculate one second derivative value
                    if is_uniform
                        dB = (B0(2) - B0(1)) * 1000; % [mT]
                        second_derivative(:) = (transition_freq(3) - 2*transition_freq(2) + transition_freq(1)) / (dB^2);
                    else
                        % Non-uniform 3-point second derivative
                        dB1 = (B0(2) - B0(1)) * 1000; % [mT]
                        dB2 = (B0(3) - B0(2)) * 1000; % [mT]
                        deriv1 = (transition_freq(2) - transition_freq(1)) / dB1;
                        deriv2 = (transition_freq(3) - transition_freq(2)) / dB2;
                        second_derivative(:) = 2 * (deriv2 - deriv1) / (dB1 + dB2);
                    end
                else
                    % Not enough points for second derivative
                    second_derivative = NaN(1, num_fields);
                end

                d2T{t_idx}(trans_idx, :) = second_derivative;
            end
        end
    end
end

% Display derivative statistics with state labels
fprintf('\nDerivative Statistics:\n');
fprintf('======================\n');
for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(dT{t_idx})
        deriv_data = dT{t_idx};
        fprintf('\nTransition order %d:\n', trans_order);
        fprintf('--------------------\n');

        for trans_idx = 1:size(deriv_data, 1)
            iState = mode_info{t_idx}.iState(trans_idx);
            fState = mode_info{t_idx}.fStates(trans_idx);

            % First derivative stats
            valid_deriv = deriv_data(trans_idx, ~isnan(deriv_data(trans_idx, :)));
            if ~isempty(valid_deriv)
                max_deriv = max(abs(valid_deriv));
                mean_deriv = mean(valid_deriv);
                std_deriv = std(valid_deriv);

                fprintf('  |%d> → |%d>:\n', iState, fState);
                fprintf('    1st derivative: Max |d/dB| = %.3e GHz/mT, Mean = %.3e ± %.3e GHz/mT\n', ...
                    max_deriv, mean_deriv, std_deriv);

                % Second derivative stats if calculated
                if Options.d2T_dB2 && ~isempty(d2T{t_idx})
                    second_deriv_data = d2T{t_idx}(trans_idx, :);
                    valid_second = second_deriv_data(~isnan(second_deriv_data));
                    if ~isempty(valid_second)
                        max_second = max(abs(valid_second));
                        mean_second = mean(valid_second);
                        std_second = std(valid_second);
                        fprintf('    2nd derivative: Max |d²/dB²| = %.3e GHz/mT², Mean = %.3e ± %.3e GHz/mT²\n', ...
                            max_second, mean_second, std_second);
                    end
                end
            end
        end
    end
end
fprintf('======================\n\n');

% Plot derivatives with state labels
if Options.dT_dB && Options.d2T_dB2
    % Plot both the first and second derivatives
    figure;
    subplot(2,1,1)
    hold on
    title('First Derivatives of Transitions Modes')
elseif Options.dT_dB && ~Options.d2T_dB2
    % Plot only the first derivatives
    figure;
    hold on
    title('First Derivatives of Transition Modes')
end

B_plot = 1000 * vecnorm(Bfield,2,1) .* sign(sum(Bfield,1)); % [mT]

% Extend line styles and colors if needed
num_transitions = length(Options.ndE);
if num_transitions > length(Options.line_styles)
    Options.line_styles = repmat(Options.line_styles, 1, ceil(num_transitions/length(Options.line_styles)));
end
if num_transitions > length(Options.colors)
    Options.colors = repmat(Options.colors, 1, ceil(num_transitions/length(Options.colors)));
end

lgd_labl = {};
pHandl = [];

for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(dT{t_idx})
        data_to_plot = dT{t_idx};
        num_trans = size(data_to_plot, 1);

        for trans_idx = 1:num_trans
            % Use different shades/styles for each transition
            color_shade = Options.colors{mod(t_idx-1, length(Options.colors))+1};
            if ischar(color_shade)
                switch color_shade
                    case 'b', base_color = [0 0 1];
                    case 'r', base_color = [1 0 0];
                    case 'g', base_color = [0 1 0];
                    case 'm', base_color = [1 0 1];
                    case 'c', base_color = [0 1 1];
                    case 'k', base_color = [0 0 0];
                    case 'y', base_color = [1 1 0];
                    otherwise, base_color = [0 0 1];
                end
                brightness_factor = 0.7 + 0.3 * (trans_idx-1)/(max(num_trans-1,1));
                plot_color = base_color * brightness_factor;
                plot_color(plot_color > 1) = 1;
            else
                plot_color = color_shade;
            end

            linStyl = Options.line_styles{mod(trans_idx-1, length(Options.line_styles))+1};

            p = plot(B_plot, data_to_plot(trans_idx, :), linStyl, ...
                'Color', plot_color, 'LineWidth', 1.5);

            iState = mode_info{t_idx}.iState(trans_idx);
            fState = mode_info{t_idx}.fStates(trans_idx);
            label = sprintf('d/dB |%d⟩→|%d⟩', iState, fState);

            pHandl = [pHandl, p];
            lgd_labl{end+1} = label;
        end
    end
end

xlabel('Magnetic Field (mT)')
ylabel('dE/dB (GHz/mT)')
grid on
if ~isempty(pHandl)
    legend(pHandl, lgd_labl, 'Location', 'bestoutside', 'FontSize', 9)
end

% Plot second derivatives if calculated
if Options.d2T_dB2 && ~isempty(d2T)
    subplot(2,1,2)
    hold on

    lgd_labl = {};
    pHandl = [];
    for t_idx = 1:length(Options.ndE)
        trans_order = Options.ndE(t_idx);
        if trans_order < En_max && ~isempty(d2T{t_idx})
            data_to_plot = d2T{t_idx};
            num_trans = size(data_to_plot, 1);

            for trans_idx = 1:num_trans
                % Use different shades/styles for each transition
                color_shade = Options.colors{mod(t_idx-1, length(Options.colors))+1};
                if ischar(color_shade)
                    switch color_shade
                        case 'b', base_color = [0 0 1];
                        case 'r', base_color = [1 0 0];
                        case 'g', base_color = [0 1 0];
                        case 'm', base_color = [1 0 1];
                        case 'c', base_color = [0 1 1];
                        case 'k', base_color = [0 0 0];
                        case 'y', base_color = [1 1 0];
                        otherwise, base_color = [0 0 1];
                    end
                    brightness_factor = 0.7 + 0.3 * (trans_idx-1)/(max(num_trans-1,1));
                    plot_color = base_color * brightness_factor;
                    plot_color(plot_color > 1) = 1;
                else
                    plot_color = color_shade;
                end

                linStyl = Options.line_styles{mod(trans_idx-1, length(Options.line_styles))+1};

                p = plot(B_plot, data_to_plot(trans_idx, :), linStyl, ...
                    'Color', plot_color, 'LineWidth', 1.5);

                iState = mode_info{t_idx}.iState(trans_idx);
                fState = mode_info{t_idx}.fStates(trans_idx);
                label = sprintf('d²/dB² |%d⟩ → |%d⟩', iState, fState);

                pHandl = [pHandl, p];
                lgd_labl{end+1} = label;
            end
        end
    end

    xlabel('Magnetic Field (mT)')
    ylabel('d²E/dB² (GHz/mT²)')
    title('Second Derivatives of Transitions (with State Labels)')
    grid on
    if ~isempty(pHandl)
        legend(pHandl, lgd_labl, 'Location', 'bestoutside', 'FontSize', 9)
    end
end

% ZEFOZ Filtering Section - Filter transitions by derivative criteria
if isfield(Options, 'zefoz') && length(Options.zefoz) >= 3
    fprintf('\n=== ZEFOZ Filtering ===\n');
    fprintf('Field window: %.3f - %.3f mT\n', Options.zefoz(1)*1000, Options.zefoz(2)*1000);
    fprintf('Derivative threshold: ±%.2e GHz/mT\n', Options.zefoz(3));

    % Initialize storage for filtered transitions
    zefoz_modes = struct();
    zefoz_modes.mode_order = [];
    zefoz_modes.mode_idx = [];
    zefoz_modes.iState = [];
    zefoz_modes.fState = [];
    zefoz_modes.freq = {};  % Full frequency array for each filtered transition
    zefoz_modes.dT_dB = {};  % Full derivative array for each filtered transition
    zefoz_modes.d2T_dB2 = {};  % Full second derivative array for each filtered transition

    % Define field window
    B_range = sort(Options.zefoz(1:2));
    B_low = B_range(1); % [T]
    B_hi = B_range(2); % [T]
    zefoz_criteria = Options.zefoz(3); % [GHz/mT]

    % Find indices within field window
    B0 = vecnorm(Bfield,2,1) .* sign(sum(Bfield,1));
    field_mask = B0 >= B_low & B0 <= B_hi;
    field_indices = find(field_mask);

    if isempty(field_indices)
        warning('No field points found within the specified window');
    else
        fprintf('Found %d field points within window\n', length(field_indices));

        % Check each transition
        for t_idx = 1:length(Options.ndE)
            trans_order = Options.ndE(t_idx);
            if trans_order < En_max && ~isempty(dT{t_idx})
                deriv_data = dT{t_idx};
                freq_data = deltaEn{t_idx};

                % Get second derivative data if available
                if Options.d2T_dB2 && ~isempty(d2T{t_idx})
                    second_deriv_data = d2T{t_idx};
                else
                    second_deriv_data = [];
                end

                for trans_idx = 1:size(deriv_data, 1)
                    % Get derivatives within field window
                    deriv_in_window = deriv_data(trans_idx, field_indices);

                    % Check if ALL derivatives within window are below threshold
                    if all(abs(deriv_in_window) <= abs(zefoz_criteria))
                        % This transition meets the criteria
                        zefoz_modes.mode_order(end+1) = trans_order;
                        zefoz_modes.mode_idx(end+1) = trans_idx;
                        zefoz_modes.iState(end+1) = mode_info{t_idx}.iState(trans_idx);
                        zefoz_modes.fState(end+1) = mode_info{t_idx}.fStates(trans_idx);

                        % Store the FULL frequency and derivative arrays for this transition
                        zefoz_modes.freq{end+1} = freq_data(trans_idx, :); % Full array [GHz]
                        zefoz_modes.dT_dB{end+1} = deriv_data(trans_idx, :); % Full array [GHz/mT]

                        % Calculate second derivative if not already computed
                        if ~isempty(second_deriv_data)
                            % Use existing second derivative
                            zefoz_modes.d2T_dB2{end+1} = second_deriv_data(trans_idx, :);
                        else
                            % Calculate second derivative for this specific transition
                            transition_freq = freq_data(trans_idx, :);
                            num_fields = length(transition_freq);
                            second_derivative = zeros(1, num_fields);

                            if num_fields > 3
                                if is_uniform
                                    dB = B0(2) - B0(1); % [T]
                                    dB2 = (dB * 1000)^2; % Convert to (mT)^2

                                    % Forward difference for first point
                                    second_derivative(1) = (2*transition_freq(1) - 5*transition_freq(2) + 4*transition_freq(3) - transition_freq(4)) / dB2;

                                    % Central difference for interior points
                                    for ii = 2:num_fields-1
                                        second_derivative(ii) = (transition_freq(ii+1) - 2*transition_freq(ii) + transition_freq(ii-1)) / dB2;
                                    end

                                    % Backward difference for last point
                                    second_derivative(end) = (-transition_freq(end-3) + 4*transition_freq(end-2) - 5*transition_freq(end-1) + 2*transition_freq(end)) / dB2;
                                else
                                    % Non-uniform grid handling
                                    for ii = 2:num_fields-1
                                        dB_back = (B0(ii) - B0(ii-1)) * 1000; % [mT]
                                        dB_forw = (B0(ii+1) - B0(ii)) * 1000; % [mT]

                                        deriv_back = (transition_freq(ii) - transition_freq(ii-1)) / dB_back;
                                        deriv_forw = (transition_freq(ii+1) - transition_freq(ii)) / dB_forw;

                                        second_derivative(ii) = 2 * (deriv_forw - deriv_back) / (dB_back + dB_forw);
                                    end
                                    % Handle boundaries
                                    second_derivative(1) = second_derivative(2);
                                    second_derivative(end) = second_derivative(end-1);
                                end
                            else
                                second_derivative = NaN(1, num_fields);
                            end

                            zefoz_modes.d2T_dB2{end+1} = second_derivative;
                        end

                        % Calculate statistics for this transition (using local variables only)
                        freq_in_window = freq_data(trans_idx, field_indices);
                        mean_deriv = mean(abs(deriv_in_window));
                        max_deriv = max(abs(deriv_in_window));
                        mean_freq = mean(freq_in_window);

                        % Calculate second derivative statistics
                        second_deriv_window = zefoz_modes.d2T_dB2{end}(field_indices);
                        valid_second = second_deriv_window(~isnan(second_deriv_window));
                        if ~isempty(valid_second)
                            mean_d2 = mean(abs(valid_second));
                            max_d2 = max(abs(valid_second));
                        else
                            mean_d2 = NaN;
                            max_d2 = NaN;
                        end

                        fprintf('  Found: |%d⟩ → |%d⟩ (Order %d, Trans %d)\n', ...
                            zefoz_modes.iState(end), ...
                            zefoz_modes.fState(end), ...
                            trans_order, trans_idx);
                        fprintf('    Mean |d/dB| = %.3e GHz/mT, Max |d/dB| = %.3e GHz/mT\n', ...
                            mean_deriv, max_deriv);
                        fprintf('    Mean |d²/dB²| = %.3e GHz/mT², Max |d²/dB²| = %.3e GHz/mT²\n', ...
                            mean_d2, max_d2);
                        fprintf('    Mean frequency = %.3f GHz\n', mean_freq);
                    end
                end
            end
        end

        fprintf('\nTotal filtered transitions: %d\n', length(zefoz_modes.mode_order));
    end

    % Plot filtered transitions with three subplots
    if ~isempty(zefoz_modes.mode_order)
        figure;

        % First subplot: Transition frequencies
        subplot(3,1,1)
        hold on
        box on

        % Plot the field window as shaded region
        y_lim_temp = [0 1]; % Temporary, will adjust after plotting
        patch([B_low B_hi B_hi B_low]*1000, ...
            [y_lim_temp(1) y_lim_temp(1) y_lim_temp(2) y_lim_temp(2)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

        % Plot each filtered transition's frequency
        pHandl = [];
        lgd_labl = {};
        colors = lines(length(zefoz_modes.mode_order));

        B_plot = 1000 * vecnorm(Bfield,2,1) .* sign(sum(Bfield,1)); % [mT]
        for i = 1:length(zefoz_modes.mode_order)
            % Plot full transition curve using freq_full
            plot(B_plot, zefoz_modes.freq{i}, '-', ...
                'Color', [colors(i,:) 0.3], 'LineWidth', 1);

            % Highlight the portion within the field window
            p = plot(B_plot(field_indices), ...
                zefoz_modes.freq{i}(field_indices), ...
                '-', 'Color', colors(i,:), 'LineWidth', 2);

            pHandl = [pHandl, p];
            lgd_labl{end+1} = sprintf('|%d⟩→|%d⟩', ...
                zefoz_modes.iState(i), ...
                zefoz_modes.fState(i));
        end

        % Adjust y-limits based on plotted data
        ylim_vals = ylim;
        delete(findobj(gca, 'Type', 'patch')); % Remove old patch
        patch([B_low B_hi B_hi B_low]*1000, ...
            [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

        % Replot to ensure data is on top of patch
        uistack(findobj(gca, 'Type', 'line'), 'top');

        xlabel('Magnetic Field (mT)')
        ylabel('Frequency (GHz)')
        title(sprintf('ZEFOZ Transitions (|dE/dB| ≤ %.3e GHz/mT in [%.1f, %.1f] mT)', ...
            zefoz_criteria, B_low*1000, B_hi*1000))
        grid on
        legend(pHandl, lgd_labl, 'Location', 'best', 'FontSize', 9)

        % Second subplot: First derivatives
        subplot(3,1,2)
        hold on
        box on

        % Plot the field window and derivative threshold
        patch([B_low B_hi B_hi B_low]*1000, ...
            [-abs(zefoz_criteria) -abs(zefoz_criteria) abs(zefoz_criteria) abs(zefoz_criteria)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'k', 'FaceAlpha', 0.3, 'LineStyle', '--');
    
        B_plot = 1000 * vecnorm(Bfield,2,1) .* sign(sum(Bfield,1)); % [mT]
        for i = 1:length(zefoz_modes.mode_order)
            % Plot full derivative curve using dT_dB_full
            plot(B_plot, zefoz_modes.dT_dB{i}, '-', ...
                'Color', [colors(i,:) 0.3], 'LineWidth', 1);

            % Highlight the portion within the field window
            plot(B_plot(field_indices), ...
                zefoz_modes.dT_dB{i}(field_indices), ...
                '-', 'Color', colors(i,:), 'LineWidth', 2);
        end

        % Add horizontal line at y=0
        plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);

        xlabel('Magnetic Field (mT)')
        ylabel('dE/dB (GHz/mT)')
        title('First Derivatives')
        grid on

        % Third subplot: Second derivatives
        subplot(3,1,3)
        hold on
        box on

        % Plot the field window as shaded region
        y_lim_temp = ylim;
        patch([B_low B_hi B_hi B_low]*1000, ...
            [y_lim_temp(1) y_lim_temp(1) y_lim_temp(2) y_lim_temp(2)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

        B_plot = 1000 * vecnorm(Bfield,2,1) .* sign(sum(Bfield,1)); % [mT]
        for i = 1:length(zefoz_modes.mode_order)
            % Plot full second derivative curve
            plot(B_plot(3,:)*1000, zefoz_modes.d2T_dB2{i}, '-', ...
                'Color', [colors(i,:) 0.3], 'LineWidth', 1);

            % Highlight the portion within the field window
            plot(B_plot(field_indices), ...
                zefoz_modes.d2T_dB2{i}(field_indices), ...
                '-', 'Color', colors(i,:), 'LineWidth', 2);
        end

        % Adjust y-limits and redraw patch
        ylim_vals = ylim;
        delete(findobj(gca, 'Type', 'patch'));
        patch([B_low B_hi B_hi B_low]*1000, ...
            [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        uistack(findobj(gca, 'Type', 'line'), 'top');

        % Add horizontal line at y=0
        plot(xlim, [0 0], 'k-', 'LineWidth', 0.5);

        xlabel('Magnetic Field (mT)')
        ylabel('d²E/dB² (GHz/mT²)')
        title('Second Derivatives')
        grid on

        sgtitle(sprintf('ZEFOZ Analysis: %d Transitions Meeting Criteria', ...
            length(zefoz_modes.mode_order)), 'FontSize', 12, 'FontWeight', 'bold')
    else
        fprintf('\nNo transitions found meeting ZEFOZ criteria.\n');
        zefoz_modes = [];
    end

    % Store filtered transitions in options for output
    Options.zefoz_modes = zefoz_modes;

else
    if isfield(Options, 'zefoz')
        warning('Options.zefoz must have at least 3 elements: [B_low, B_hi, derivative_threshold]');
    end
    zefoz_modes = []; % Return empty if no ZEFOZ analysis
end
end
