%% Main Script
clearvars

% Setup options 
Options.CEF = false; % use crystal field parameters
Options.nEn = 16; % number of eigenstates to include
Options.ndE = [1:16]; % specify which transition orders to calculate and plot
% Options.ndE = [4, 6, 8, 10, 11, 12, 13, 14, 15]; % specify which transition orders to calculate and plot
Options.dT_dB = false; % First derivatives of transitions
Options.d2T_dB2 = false; % Second derivatives of transitions
Options.chi = false; % calculate magnetic susceptibility
Options.freq = linspace(0, 4, 1001); % frequency range for susceptibility [GHz]
Options.temp = 0.1; % temperature for susceptibility calculation [K]
Options.gamma = 1e-5; % spin linewidth [meV]

Options.line_styles = {'-', '--', '-.', ':'}; % line styles for different transitions
Options.colors = {'b', 'r', 'g', 'm', 'c', 'k', 'y'}; % colors for different transitions

B0 = linspace(-0.2, 0.2, 801); % [mT]
theta = 0; % [degree] deviation angle from c-axis
phi = 0; % [degree] rotation angle within ab-plane, relative to a-axis
Bfield = [B0*sin(pi*theta/180)*cos(pi*phi/180); B0*sin(pi*theta/180)*sin(pi*phi/180); B0*cos(pi*theta/180)];

% Initialize system parameters and operators
[const, params, ham_cf] = initialization(Options);
% En_max = (2*params.J+1)*(2*params.I+1); % maximal number of states

% Calculate eigenvalues and eigenvectors
[eigenE, eigenW, dEn] = eigenEnergy(Bfield, Options, const, params, ham_cf);

% Calculate transitions and derivatives
deltaE(Options, params.J, params.I, dEn, Bfield);

% Calculate susceptibility if requested
if Options.chi
    results = MF_chi(Options, const, params.J, params.I, eigenE, eigenW, Bfield);
end

fprintf('All calculations completed.\n');

function [const, params, ham_cf] = initialization(Options)
% define constants
const.hbar = 1.05457E-34; % Reduced Planck constant [J.s]
const.muB = 9.274e-24; % [J/T]
const.kB = 1.3806e-23; % [J/K]
const.muN = 5.05078e-27; % [J/T]
const.mu0 = 4e-7 * pi; % [H/m]
const.kB_meV = 8.61733e-2; % Boltzmann constant [meV/K]
const.J2meV = 6.24151e+21; % [mev/J]
const.Gh2mV = const.hbar * 2*pi * 10^9 * const.J2meV; % [meV/GHz]

I = 7/2;
if Options.CEF == true
    % Option 1: use CEF
    J = 15/2; % Er
    Bxx = [567.02 -945.07 1008.10 0 -22.16  936.29  0.854] * 0.123983; % [meV] J. Chem. Phys. 53 (9), 1970
    % B = [548.81 -942.41 1004.66 0 -17.12  946.84  1.343] * 0.123983; % [meV] J. Chem. Phys. 53 (9), 1970
    % B = [238 -90 -852 0 -0.6 -396 -75] * 0.123983; % [meV] J. Chem. Phys. 55, 5 (1971)
    % B = [640 -809 860 0 12  444  161] * 0.123983; % [meV] J. Chem. Phys. 54 314, (1971)
    % B = [732 -851 963 0 -9  555  105] * 0.123983; % [meV] J. Chem. Phys. 55, 2538 (1971)
    % B = [292 -832 834 0 -62  445  92] * 0.123983; % [meV] J. Chem. Phys. 54 314, (1971)
    % B = [308 -904 906 0 -77 558  115] * 0.123983; % [meV] J. Chem. Phys. 55, 2538 (1971)
    ham_cf = cf(J, Bxx, 0);
    ham_cf = kron(ham_cf, eye(2*I+1)); % Expand the crystal field space to include nuclear moments' degrees of freedom
else
    % Option 2: use effective doublet hamiltonian
    J = 1/2; % effective Kramer's spin
    ham_cf = 0;
end
% hyperfine parameters
A = [-8.71436227E2 -8.71436227E2 1.10543961E2] / 1000 * const.Gh2mV; % [meV]
% A = -[873 873 130] / 1000 * const.Gh2mV; % [meV]
% A = [-8.87920245E2 -8.87920245E2 1.67888823E2] / 1000 * const.Gh2mV; % [meV]

% electronic g-factor
gE = [8.38 8.38 1.26234449];
% gE = [8.38 8.38 1.247];

gN = [-1.36503716E-1 -1.36503716E-1 -1.36503716E-1]; % nuclear g-factor
Q = [-2.55041373 -2.55041373 -5.28392734] / 1000 * const.Gh2mV; % nuclear quadrupler interaction strength

params.J = J;
params.I = I;
params.A = A;
params.gE = gE;
params.gN = gN;
params.Q = Q;
end

function [eigenE, eigenW, dEn] = eigenEnergy(Bfield, Options, const, params, ham_cf)
J = params.J;
I = params.I;
A = params.A;
Q = params.Q;
gE = params.gE;
gN = params.gN;

En_max = (2*J+1)*(2*I+1); % maximal number of states
if Options.nEn > En_max; Options.nEn = En_max; end

[~,~,~,~,~,~,Jxh,Jyh,Jzh,Ixh,Iyh,Izh] = spin_operators(J, I);

% Initialize arrays
eigenE = zeros(Options.nEn, size(Bfield,2));
eigenW = zeros(Options.nEn, Options.nEn, size(Bfield,2));
dEn = zeros(Options.nEn, Options.nEn, size(Bfield,2));
for ii = 1:size(Bfield,2)
    ham_E = const.muB * const.J2meV *...
        (gE(1) * Jxh * Bfield(1,ii) + gE(2) * Jyh * Bfield(2,ii) + gE(3) * Jzh * Bfield(3,ii));
    % ham_N = A(1)*Ixh*Jxh + A(2)*Iyh*Jyh + A(3)*Izh*Jzh;
    ham_N = A(1)*Ixh*Jxh + A(2)*Iyh*Jyh + A(3)*Izh*Jzh +...
        Q(1)*Ixh*Ixh + Q(2)*Iyh*Iyh + Q(3)*Izh*Izh + const.muN *...
        (gN(1) * Ixh * Bfield(1,ii) + gN(2) * Iyh * Bfield(2,ii) + gN(3) * Izh * Bfield(3,ii));

    ham = ham_cf + ham_E + ham_N;

    [wv, ee] = eig(ham);
    ee = real(diag(ee)); % Take only the real part of the eigen-energy to form a diagonal matrix
    [ee, n] = sort(ee); % sort the energy from lowest to the highest
    wv = wv(:,n);
    ee = ee - min(ee); % Normalize the energy amplitude to the lowest eigen-energy
    eigenE(:,ii) = ee(1:Options.nEn);
    eigenW(:,:,ii) = wv(1:Options.nEn,1:Options.nEn); % sort the eigen-vectors in its basis accordingly

    % calculate all the possible inter-state transitions
    en = ee(1:Options.nEn) / const.Gh2mV; % [GHz]
    [Ex, Ey] = meshgrid(en,en);
    dEn(:,:,ii) = Ex - Ey; % all possible transitions among eigenstates
end

end

function deltaE(Options, J, I, dEn, Bfield)
En_max = (2*J+1)*(2*I+1); % maximal number of states
if Options.nEn > En_max; Options.nEn = En_max; end
% Calculate transitions for specified orders
deltaEn = cell(length(Options.ndE), 1);

for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < Options.nEn
        % Initialize array for this transition order
        num_transitions = Options.nEn - trans_order;
        deltaEn{t_idx} = zeros(num_transitions, size(dEn,3));

        % Calculate transitions for each magnetic field point
        for ii = 1:size(dEn,3)
            deltaEn{t_idx}(:, ii) = diag(squeeze(dEn(:,:,ii)), trans_order);
        end
    end
end

fprintf('Calculated transitions for orders: %s\n', mat2str(Options.ndE));
fprintf('Total number of eigenstates: %d\n', Options.nEn);

% Plot basic transitions
plot_basic_transitions(Options, deltaEn, Bfield, En_max);

% Calculate derivatives if requested
if Options.dT_dB
    dT_dB(Options, deltaEn, Bfield, En_max);
end
end

function plot_basic_transitions(Options, deltaEn, Bfield, En_max)
% Plotting
EigenModes = figure;
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

Bz = Bfield(3,:); % [T]

for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(deltaEn{t_idx})
        % Get the data for this transition order
        data_to_plot = deltaEn{t_idx};

        if ~isempty(data_to_plot)
            p = plot(Bz*1000, data_to_plot, Options.line_styles{t_idx}, ...
                'Color', Options.colors{t_idx}, 'LineWidth', 1.5);
            plot_handles = [plot_handles, p(1)];
            legend_labels{end+1} = sprintf('|n\\rangle \\rightarrow |n+%d\\rangle', trans_order);
        end
    end
end

% Customize plot
xlabel('Magnetic Field (mT) $\|$ c','interpreter','latex')
ylabel('Energy gap (GHz)')
set(gca,'fontsize',14)

% Add legend
if ~isempty(plot_handles)
    legend(plot_handles, legend_labels, 'Location', 'best')
end

grid on
title('Eigenstate Transitions vs Magnetic Field')
end

function dT_dB(Options, deltaEn, Bfield, En_max)
% Calculate first and second derivatives of transitions
fprintf('\nCalculating transition derivatives...\n');

Bz = Bfield(3,:); % [T]
transition_derivatives = cell(length(Options.ndE), 1);
transition_second_derivatives = cell(length(Options.ndE), 1);

for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(deltaEn{t_idx})
        data_to_process = deltaEn{t_idx};
        num_transitions = size(data_to_process, 1);
        num_fields = size(data_to_process, 2);

        % Initialize derivative arrays
        transition_derivatives{t_idx} = zeros(num_transitions, num_fields);
        if Options.d2T_dB2
            transition_second_derivatives{t_idx} = zeros(num_transitions, num_fields);
        end

        % Calculate derivatives for each transition
        for trans_idx = 1:num_transitions
            transition_freq = data_to_process(trans_idx, :);

            % First derivative calculation
            derivative = zeros(1, num_fields);

            % Forward difference for first point
            if num_fields > 1
                dB = Bz(2) - Bz(1); % [T]
                derivative(1) = (transition_freq(2) - transition_freq(1)) / (dB * 1000); % [GHz/mT]
            end

            % Central difference for middle points
            for ii = 2:num_fields-1
                dB_forward = Bz(ii+1) - Bz(ii); % [T]
                dB_backward = Bz(ii) - Bz(ii-1); % [T]
                dB_avg = (dB_forward + dB_backward) / 2; % [T]

                derivative(ii) = (transition_freq(ii+1) - transition_freq(ii-1)) / (2 * dB_avg * 1000); % [GHz/mT]
            end

            % Backward difference for last point
            if num_fields > 1
                dB = Bz(end) - Bz(end-1); % [T]
                derivative(end) = (transition_freq(end) - transition_freq(end-1)) / (dB * 1000); % [GHz/mT]
            end

            transition_derivatives{t_idx}(trans_idx, :) = derivative;

            % Second derivative calculation
            if Options.d2T_dB2 && num_fields > 2
                second_derivative = zeros(1, num_fields);

                % Forward difference for first point
                dB = Bz(2) - Bz(1); % [T]
                second_derivative(1) = (transition_freq(3) - 2*transition_freq(2) + transition_freq(1)) / ((dB * 1000)^2);

                % Central difference for middle points
                for ii = 2:num_fields-1
                    dB_avg = (Bz(ii+1) - Bz(ii-1)) / 2; % [T]
                    second_derivative(ii) = (transition_freq(ii+1) - 2*transition_freq(ii) + transition_freq(ii-1)) / ((dB_avg * 1000)^2);
                end

                % Backward difference for last point
                dB = Bz(end) - Bz(end-1); % [T]
                second_derivative(end) = (transition_freq(end) - 2*transition_freq(end-1) + transition_freq(end-2)) / ((dB * 1000)^2);

                transition_second_derivatives{t_idx}(trans_idx, :) = second_derivative;
            end
        end
    end
end

% Display derivative statistics
fprintf('\nDerivative Statistics:\n');
for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(transition_derivatives{t_idx})
        deriv_data = transition_derivatives{t_idx};
        fprintf('Transition order %d:\n', trans_order);
        for trans_idx = 1:size(deriv_data, 1)
            max_deriv = max(abs(deriv_data(trans_idx, :)));
            mean_deriv = mean(deriv_data(trans_idx, :));
            fprintf('  Transition %d: Max |d/dB| = %.4f GHz/mT, Mean d/dB = %.4f GHz/mT\n', ...
                trans_idx, max_deriv, mean_deriv);
        end
    end
end

% Plot derivatives if requested
if Options.dT_dB && Options.d2T_dB2
    % Plot both the first and second derivatives
    figure;
    subplot(2,1,1)
    hold on
elseif Options.dT_dB && ~Options.d2T_dB2
    % Plot only the first derivatives
    figure;
    hold on
end
Bz = Bfield(3,:) * 1000; % [mT]

% Extend line styles and colors if needed
num_transitions = length(Options.ndE);
if num_transitions > length(Options.line_styles)
    Options.line_styles = repmat(Options.line_styles, 1, ceil(num_transitions/length(Options.line_styles)));
end
if num_transitions > length(Options.colors)
    Options.colors = repmat(Options.colors, 1, ceil(num_transitions/length(Options.colors)));
end

for t_idx = 1:length(Options.ndE)
    trans_order = Options.ndE(t_idx);
    if trans_order < En_max && ~isempty(transition_derivatives{t_idx})
        data_to_plot = transition_derivatives{t_idx};
        plot(Bz, data_to_plot, Options.line_styles{t_idx}, ...
            'Color', Options.colors{t_idx}, 'LineWidth', 1.5);
    end
end

xlabel('Magnetic Field (mT)')
ylabel('dE/dB (GHz/mT)')
title('First Derivatives of Transitions')
grid on

% Plot second derivatives if calculated
if Options.d2T_dB2 && ~isempty(transition_second_derivatives)
    subplot(2,1,2)
    hold on

    for t_idx = 1:length(Options.ndE)
        trans_order = Options.ndE(t_idx);
        if trans_order < En_max && ~isempty(transition_second_derivatives{t_idx})
            data_to_plot = transition_second_derivatives{t_idx};
            plot(Bz, data_to_plot, Options.line_styles{t_idx}, ...
                'Color', Options.colors{t_idx}, 'LineWidth', 1.5);
        end
    end

    xlabel('Magnetic Field (mT)')
    ylabel('d²E/dB² (GHz/mT²)')
    title('Second Derivatives of Transitions')
    grid on
end
end