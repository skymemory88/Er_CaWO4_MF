% ndE = Options.ndE;
ndE = 3; % transition order
plot_handles = [];
% figure
hold on

% Calculate transitions for all orders first
for t_idx = 1:length(ndE)
    trans_order = ndE(t_idx);
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

legend_labels = {};

% Calculate maximum number of transitions in any group to determine color palette size
max_transitions = 0;
for t_idx = 1:length(ndE)
    trans_order = ndE(t_idx);
    if trans_order < Options.nEn
        num_transitions = Options.nEn - trans_order;
        max_transitions = max(max_transitions, num_transitions);
    end
end

colors = lines(max_transitions); % Colors for transitions within each group
line_styles = {'-', '--', '-.', ':', '-', '--', '-.'}; % Line styles for different groups

Bz = Bfield(3,:); % [T]

% Plot each individual transition
for t_idx = 1:length(ndE)
    trans_order = ndE(t_idx);
    if trans_order < Options.nEn && ~isempty(deltaEn{t_idx})
        
        % Get the data for this transition order
        data_to_plot = deltaEn{t_idx};
        num_transitions = size(data_to_plot, 1);
        
        % Get line style for this transition order group
        style_idx = mod(t_idx-1, length(line_styles)) + 1;
        current_line_style = line_styles{style_idx};
        
        % Plot each individual transition |n> -> |n+t>
        for n = 1:num_transitions
            initial_state = n - 1;  % Convert to 0-based indexing for state labels
            final_state = initial_state + trans_order;
            
            % Use different colors within each transition order group
            color_idx = mod(n-1, size(colors,1)) + 1;
            
            p = plot(Bz, data_to_plot(n, :), current_line_style, ...
                'Color', colors(color_idx,:), 'LineWidth', 1.5);
            
            plot_handles = [plot_handles, p];
            legend_labels{end+1} = sprintf('|%d\\rangle \\rightarrow |%d\\rangle', ...
                initial_state, final_state);
        end
    end
end

% Add the legend
legend(plot_handles, legend_labels, 'Location', 'best', 'Interpreter', 'tex');

% Optional: Group legend entries by transition order for better organization
% Uncomment the following section if you prefer grouped legend entries
% % Alternative: Group by transition order with separators
% legend_labels_grouped = {};
% plot_handles_grouped = [];
% 
% for t_idx = 1:length(ndE)
%     trans_order = ndE(t_idx);
%     if trans_order < Options.nEn && ~isempty(deltaEn{t_idx})
%         
%         % Add a group header (you might need to create a dummy invisible plot for this)
%         if t_idx > 1
%             % Add some spacing in legend
%             legend_labels_grouped{end+1} = '';
%         end
%         
%         data_to_plot = deltaEn{t_idx};
%         num_transitions = size(data_to_plot, 1);
%         
%         % Get line style for this transition order group
%         style_idx = mod(t_idx-1, length(line_styles)) + 1;
%         current_line_style = line_styles{style_idx};
%         
%         for n = 1:num_transitions
%             initial_state = n - 1;
%             final_state = initial_state + trans_order;
%             
%             % Use different colors within each transition order group
%             color_idx = mod(n-1, size(colors,1)) + 1;
%             
%             p = plot(Bz, data_to_plot(n, :), current_line_style, ...
%                 'Color', colors(color_idx,:), 'LineWidth', 1.5);
%             
%             plot_handles_grouped = [plot_handles_grouped, p];
%             legend_labels_grouped{end+1} = sprintf('|%d\\rangle \\rightarrow |%d\\rangle', ...
%                 initial_state, final_state);
%         end
%     end
% end