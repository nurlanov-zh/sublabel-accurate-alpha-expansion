function [ result_sublabel, time_elapsed_sublabel, energy_raw_sublabel, ...
    valueD, valueS, count_non_convex] = sublabel_local_QL(label_inds_gco, ...
            label_space_double, data_cost, smooth_cost, neighbor_rows, ...
            neighbor_cols)
%   Initialize by GCO, solve sublabel accurate convex problem by quadratic
%   local interpolation of data term and linear local interpolation of
%   smoothness term

    % consider using gurobi solver or mosek solver
    
    fprintf('\nStarting sublabel refinement...\n');
    yalmip('clear');
    
    [L, N] = size(data_cost);

    %% Data term
    tic;
    %%%%%%%%%%%%%% Neighborhood in 1D label space %%%%%%%%%%%%%%%
    label_neighbor_size = 3;
    bound_for_adaptive = 3;         % It is possible to approximate larger neighborhood
    
    [idx_matrix, quadratic, count_non_convex...
        ] = local_quadratic_data_approx(label_inds_gco, label_space_double, ...
        data_cost, label_neighbor_size, bound_for_adaptive);

    labels_matrix = zeros(N, label_neighbor_size);
    for j = 1:label_neighbor_size
        for i = 1:N
            labels_matrix(i, j) = label_space_double(idx_matrix(i, j));
        end
    end
    var_left = labels_matrix(:, 1);
    var_right = labels_matrix(:, end);

    % variable X represents the actual label of vertex
    X = sdpvar(N, 1, 'full');
    C = [X(:) >= var_left(:), X(:) <= var_right(:)];
    D = sum(X(:).*X(:).*quadratic(:, 1) + X(:).*quadratic(:, 2) + quadratic(:, 3));

    time_first_constraint = toc;
    fprintf('Time: prepare Data Term = %f sec\n', time_first_constraint);

    %% Regularization term
    % Quadratic regularization term without relaxation
    % E = E + sum(sum((X(neighbor_rows, :))' * X(neighbor_cols, :) .* smooth_cost));
    
    tic;
    
    [~, num_edges] = size(neighbor_rows);  
    
    x_left_idx = idx_matrix(neighbor_rows, :);
    x_right_idx = idx_matrix(neighbor_cols, :);
    x_left_label = labels_matrix(neighbor_rows, :);
    x_right_label = labels_matrix(neighbor_cols, :);
    approx_smooth_cost_coefs = zeros(num_edges, 1);
    approx_smooth_cost_residual = 0;
 
    for edge_ij = 1:num_edges
        [grid_m_idx, grid_n_idx] = ndgrid(x_left_idx(edge_ij, :), x_right_idx(edge_ij, :));
        pairs_idx = [grid_m_idx(:),grid_n_idx(:)];
        smooth_cost_vector = smooth_cost(sub2ind([L, L], pairs_idx(:, 1), pairs_idx(:, 2)));
        
        [grid_m_label, grid_n_label] = ndgrid(x_left_label(edge_ij, :), x_right_label(edge_ij, :));
        pairs_label = [grid_m_label(:),grid_n_label(:)];
        
        %%%%%%%% L1 norm approximation %%%%%%%%%%%%%%%
        diffs = abs(pairs_label(:, 1) - pairs_label(:, 2));
        
        %%%%%%%%% Quadratic approximation %%%%%%%%%%%%%%
        % diffs = (pairs_label(:, 1) - pairs_label(:, 2)).^2;
        
        [coefs, ~] = linsolve(diffs, smooth_cost_vector);
        
        residual = norm(diffs * coefs - smooth_cost_vector) / norm(smooth_cost_vector);
        approx_smooth_cost_residual = approx_smooth_cost_residual + residual;
        approx_smooth_cost_coefs(edge_ij, :) = coefs;
    end
    approx_smooth_cost_residual = approx_smooth_cost_residual / num_edges;
    % fprintf('approx_smooth_cost_residual = %f percent\n', approx_smooth_cost_residual);
    
    
    %%%%% L1 norm approximation via YALMIP's linear graph representations %%%%%%%%
    S = norm(approx_smooth_cost_coefs .* (X(neighbor_rows, :) - X(neighbor_cols, :)), 1);
    % S = norm((X(neighbor_rows, :) - X(neighbor_cols, :)), 1);
    
    %%%%% Quadratic approximation %%%%%%%%
    %S = sum(approx_smooth_cost_coefs .* ((X(neighbor_rows, :) - ...
    %    X(neighbor_cols, :)) .* (X(neighbor_rows, :) - X(neighbor_cols, :))));
    
    %S = norm(approx_smooth_cost_coefs .* (X(neighbor_rows, :) - X(neighbor_cols, :)), 2)^2;
    
    E = D + S;
    
    time_smooth_energy = toc;
    fprintf('Time: prepare Regularizer Term Smooth energy = %f sec\n', time_smooth_energy);
    
    %% Solve the problem
    tic;
    ops = sdpsettings('solver','gurobi,mosek,*', 'verbose',0);

    % [model,recoverymodel,diagnostic,internalmodel] = export(C,E,ops);
    % A = full(model.A);
    % dlmwrite('constraints_matrix.txt',A, ' ');

    sol = optimize(C, E, ops);
    %check(C);
    
    %% Output the solution
    if sol.problem == 0
     % Extract and display value
        solution = value(X);
        time_solve = toc;
        result_sublabel = min(1, max(0, solution));
        energy_raw_sublabel = value(E);
        valueD = value(D);
        valueS = value(S);
     
        fprintf('Time: sublabel refinement solve = %f sec\n', time_solve);
        time_elapsed_sublabel = time_first_constraint + time_smooth_energy + time_solve;
    else
        display('Hmm, something went wrong!');
        sol.info
        yalmiperror(sol.problem)
    end
end

