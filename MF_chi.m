function results = MF_chi(Options, const, J, I, eigenE, eigenW, Bfield)
    fprintf('\nCalculating magnetic susceptibility...\n');

    % Setup frequency grid
    freq_total = Options.freq; % [GHz]
    omega_grid = freq_total * const.Gh2mV; % Convert to [meV]
    
    % Setup temperature grid
    if isscalar(Options.temp)
        temperatures = ones(1, size(eigenE,2)) * Options.temp;
    else
        temperatures = Options.temp;
    end
    
    % Initialize susceptibility tensor
    chi = zeros(3, 3, length(freq_total), size(eigenE,2));
    [Jx,Jy,Jz,~,~,~,Jxh,Jyh,Jzh,Ixh,Iyh,Izh] = spin_operators(J, I);
    if I > 0 % Include hyperfine coupling
        gN = 0.1; % Simplified nuclear contribution factor
        % Hybridized operators
        JhT.x = Jxh + gN * Ixh;
        JhT.y = Jyh + gN * Iyh;
        JhT.z = Jzh + gN * Izh;
    else
        JhT.x = Jx;
        JhT.y = Jy;
        JhT.z = Jz;
    end
    
    % Calculate susceptibility for each magnetic field point
    gamma = ones(size(eigenE,1)) * Options.gamma; % [meV]
    
    for ii = 1:size(eigenE,2)
        if mod(ii, 50) == 0
            fprintf('Processing field point %d/%d\n', ii, size(eigenE,2));
        end
        
        % Get eigenvectors and eigenvalues for this field
        v = squeeze(eigenW(:,:,ii)); % eigenvectors
        en = squeeze(eigenE(:,ii)); % eigenvalues [meV]
        
        % Setup thermal population
        if temperatures(ii) > 0
            beta = 1/(const.kB_meV * temperatures(ii));
            Z = sum(exp(-beta * en));
            zn = exp(-beta * en) / Z;
            [n, np] = meshgrid(zn, zn);
            NN = n - np; % population difference
        else
            % Zero temperature: only ground state occupied
            Z = zeros(size(en));
            Z(1) = 1;
            [n, np] = meshgrid(Z, Z);
            NN = n - np;
        end
        
        % Calculate matrix elements of spin operators
        ttx = v' * JhT.x * v;
        tty = v' * JhT.y * v;
        ttz = v' * JhT.z * v;
        
        % Calculate susceptibility for each frequency
        for nf = 1:length(freq_total)
            omega = omega_grid(nf); % [meV]
            
            % Energy differences
            [ee, eep] = meshgrid(en, en);
            EE = (eep - ee - omega); % [meV]
            
            % Denominator with damping
            deno = 1 ./ (EE - 1i * gamma);
            
            % Calculate susceptibility tensor components
            chi(1,1,nf,ii) = sum(sum(ttx .* ttx' .* NN .* deno));
            chi(2,2,nf,ii) = sum(sum(tty .* tty' .* NN .* deno));
            chi(3,3,nf,ii) = sum(sum(ttz .* ttz' .* NN .* deno));
            chi(1,2,nf,ii) = sum(sum(ttx .* tty' .* NN .* deno));
            chi(1,3,nf,ii) = sum(sum(ttx .* ttz' .* NN .* deno));
            chi(2,1,nf,ii) = sum(sum(tty .* ttx' .* NN .* deno));
            chi(2,3,nf,ii) = sum(sum(tty .* ttz' .* NN .* deno));
            chi(3,1,nf,ii) = sum(sum(ttz .* ttx' .* NN .* deno));
            chi(3,2,nf,ii) = sum(sum(ttz .* tty' .* NN .* deno));
        end
    end
    
    % Store results
    results.chi_tensor = chi;
    results.freq = freq_total; % [GHz]
    results.Bfield = Bfield(3,:); % [T]
    results.temp = temperatures; % [K]
    
    B0 = vecnorm(Bfield,2,1) .* sign(sum(Bfield,1)); % [T]
    chi = results.chi_tensor;
    freq_total = results.freq;
    
    % % Plot susceptibility vs magnetic field at specific frequencies
    % figure;
    % subplot(2,2,1)
    % freq_idx = round(length(freq_total)/4); % Quarter of frequency range
    % plot(Bz*1000, real(squeeze(chi_tensor(1,1,freq_idx,:))), 'r-', 'LineWidth', 2);
    % hold on
    % plot(Bz*1000, real(squeeze(chi_tensor(2,2,freq_idx,:))), 'g-', 'LineWidth', 2);
    % plot(Bz*1000, real(squeeze(chi_tensor(3,3,freq_idx,:))), 'b-', 'LineWidth', 2);
    % xlabel('Magnetic Field (mT)')
    % ylabel('Re[\chi] (arb. units)')
    % title(sprintf('Real Susceptibility at %.2f GHz', freq_total(freq_idx)))
    % legend('\chi_{xx}', '\chi_{yy}', '\chi_{zz}')
    % grid on
    % 
    % subplot(2,2,2)
    % plot(Bz*1000, imag(squeeze(chi_tensor(1,1,freq_idx,:))), 'r--', 'LineWidth', 2);
    % hold on
    % plot(Bz*1000, imag(squeeze(chi_tensor(2,2,freq_idx,:))), 'g--', 'LineWidth', 2);
    % plot(Bz*1000, imag(squeeze(chi_tensor(3,3,freq_idx,:))), 'b--', 'LineWidth', 2);
    % xlabel('Magnetic Field (mT)')
    % ylabel('Im[\chi] (arb. units)')
    % title(sprintf('Imaginary Susceptibility at %.2f GHz', freq_total(freq_idx)))
    % legend('\chi_{xx}', '\chi_{yy}', '\chi_{zz}')
    % grid on
    % 
    % % Plot susceptibility vs frequency at specific field
    % subplot(2,2,3)
    % field_idx = round(length(Bz)/2); % Middle of field range
    % plot(freq_total, real(squeeze(chi_tensor(1,1,:,field_idx))), 'r-', 'LineWidth', 2);
    % hold on
    % plot(freq_total, real(squeeze(chi_tensor(2,2,:,field_idx))), 'g-', 'LineWidth', 2);
    % plot(freq_total, real(squeeze(chi_tensor(3,3,:,field_idx))), 'b-', 'LineWidth', 2);
    % xlabel('Frequency (GHz)')
    % ylabel('Re[\chi] (arb. units)')
    % title(sprintf('Real Susceptibility at %.1f mT', Bz(field_idx)*1000))
    % legend('\chi_{xx}', '\chi_{yy}', '\chi_{zz}')
    % grid on
    % 
    % subplot(2,2,4)
    % plot(freq_total, imag(squeeze(chi_tensor(1,1,:,field_idx))), 'r--', 'LineWidth', 2);
    % hold on
    % plot(freq_total, imag(squeeze(chi_tensor(2,2,:,field_idx))), 'g--', 'LineWidth', 2);
    % plot(freq_total, imag(squeeze(chi_tensor(3,3,:,field_idx))), 'b--', 'LineWidth', 2);
    % xlabel('Frequency (GHz)')
    % ylabel('Im[\chi] (arb. units)')
    % title(sprintf('Imaginary Susceptibility at %.1f mT', Bz(field_idx)*1000))
    % legend('\chi_{xx}', '\chi_{yy}', '\chi_{zz}')
    % grid on
    
    % 2D color plot of susceptibility
    figure;
    [B_mesh, F_mesh] = meshgrid(B0*1000, freq_total);
    
    subplot(2,2,1)
    hp1 = pcolor(B_mesh, F_mesh, real(squeeze(chi(1,1,:,:))));
    set(hp1, 'EdgeColor', 'none');
    colorbar
    xlabel('Magnetic Field (mT)')
    ylabel('Frequency (GHz)')
    title('Re[\chi_{xx}]')
    
    subplot(2,2,2)
    hp2 = pcolor(B_mesh, F_mesh, imag(squeeze(chi(1,1,:,:))));
    set(hp2, 'EdgeColor', 'none');
    colorbar
    xlabel('Magnetic Field (mT)')
    ylabel('Frequency (GHz)')
    title('Im[\chi_{xx}]')
    
    subplot(2,2,3)
    hp3 = pcolor(B_mesh, F_mesh, real(squeeze(chi(3,3,:,:))));
    set(hp3, 'EdgeColor', 'none');
    colorbar
    xlabel('Magnetic Field (mT)')
    ylabel('Frequency (GHz)')
    title('Re[\chi_{zz}]')
    
    subplot(2,2,4)
    hp4 = pcolor(B_mesh, F_mesh, imag(squeeze(chi(3,3,:,:))));
    set(hp4, 'EdgeColor', 'none');
    colorbar
    xlabel('Magnetic Field (mT)')
    ylabel('Frequency (GHz)')
    title('Im[\chi_{zz}]')
    
    colormap('jet')
end