function [ result_sublabel, time_elapsed_sublabel, energy_raw_sublabel, ...
    valueD, valueS, count_non_convex] = sublabel_local_LM(label_inds_gco, ...
            label_space_double, data_cost, smooth_cost, neighbor_rows, ...
            neighbor_cols)
%   Initialize by GCO, solve sublabel accurate convex problem by linear
%   local interpolation of data term and linearized relaxed regularization 
%   term via marginalization constraints

    count_non_convex = -1;   % since here we do not check for convexity
    
    fprintf('\nStarting sublabel refinement...\n');
    yalmip('clear');
    
    [L, N] = size(data_cost);

    %% Data term
    tic;
    %%%%%%%%%%%%%% Neighborhood in 1D label space %%%%%%%%%%%%%%%
    label_neighbor_size = 3;
    
    var_left_idx = double(min(max(1, label_inds_gco - 1), ...
                            L - label_neighbor_size + 1));
    % Create matrix of possible idxs of labels for each vertex (N, label_neighbor_size)
    idx_matrix = zeros(N, label_neighbor_size);
    for i=1:label_neighbor_size
        idx_matrix(:, i) = var_left_idx + i - 1;
    end
    
    % Create data costs matrix for shortened label neighborhood
    data_cost_short = zeros(N, label_neighbor_size);
    for i = 1:N
        for j = 1:label_neighbor_size
            data_cost_short(i, j) = data_cost(idx_matrix(i, j), i);
        end
    end
    
    labels_matrix = zeros(size(idx_matrix));
    [num_sites, num_neighbors] = size(idx_matrix);
    for j = 1:num_neighbors
        for i = 1:num_sites
            labels_matrix(i, j) = label_space_double(idx_matrix(i, j));
        end
    end
    var_left = labels_matrix(:, 1);
    var_right = labels_matrix(:, end);

    % Create variable Z for introducing marginalization
    % Z of size (N, label_neighbor_size) where the row of Z - vector of probabilities
    % P(label(x_i) = l) = Z(i, l)
    Z = sdpvar(N, label_neighbor_size, 'full');

    
    C = [Z(:) >= 0, sum(Z, 2) == ones(N, 1), ...
        sum(Z .* labels_matrix, 2) >= var_left(:), ...
        sum(Z .* labels_matrix, 2) <= var_right(:)];
    
    % C = [Z(:) >= 0, sum(Z, 2) == ones(N, 1)];
    
    D = sum(sum(Z .* data_cost_short));

    time_first_constraint = toc;
    fprintf('Time: prepare Data Term = %f sec\n', time_first_constraint);

    %% Regularization term
    % Linearized relaxed regularization term via marginalization constraints
    tic;
    [~, num_edges] = size(neighbor_rows);

    Y = sdpvar(num_edges, label_neighbor_size, label_neighbor_size, 'full');

    Z_i = Z(neighbor_rows, :);
    Z_j = Z(neighbor_cols, :);
    
    C = [C, Y(:) >= 0];
    C = [C, squeeze(sum(Y, 3)) == Z_i];
    C = [C, squeeze(sum(Y, 2)) == Z_j];
    
    time_marginal = toc;
    fprintf('Time: prepare Regularizer Term Marginal Constraint = %f sec\n', time_marginal);
        
    tic;
    %%%%%%%%%%%%%% Get smooth_cost for Y (n, label_neighbor_size, label_neighbor_size)
    edges_smooth_cost = ones(num_edges, label_neighbor_size, label_neighbor_size);
    x_i = idx_matrix(neighbor_rows, :);
    x_j = idx_matrix(neighbor_cols, :);
    for edge_ij = 1:num_edges
        [grid_m, grid_n] = ndgrid(x_i(edge_ij, :), x_j(edge_ij, :));
        pairs_of_labels = [grid_m(:),grid_n(:)];
        esq = smooth_cost(sub2ind([L, L], pairs_of_labels(:, 1), pairs_of_labels(:, 2)));
        edges_smooth_cost(edge_ij, :, :) = reshape(esq, [label_neighbor_size, label_neighbor_size]);
    end
    S = sum(sum(sum(Y .* edges_smooth_cost)));
    E = D + S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time_smooth_energy = toc;
    fprintf('Time: prepare Regularizer Term Smooth energy = %f sec\n', time_smooth_energy);

    %% Solve the problem
    tic;
     % consider using gurobi solver or mosek solver
    ops = sdpsettings('solver','gurobi,mosek,*', 'verbose',0);

    %[model,recoverymodel,diagnostic,internalmodel] = export(C,E,ops);
    %A = model.A;
    %A = full(model.A);
    %dlmwrite('constraints_matrix_L_50_N_05.txt',A, ' ');
    %uni = unimodularity_check(A);
    %disp('unimodularity: ');
    %disp(uni);

    sol = optimize(C, E, ops);
    % check(C);
    
    %% Output the solution
    if sol.problem == 0
     % Extract and display value
        valueZ = value(Z);
        time_solve = toc;
        
        % Solution as the expected value
        solution = sum(valueZ .* labels_matrix, 2);
        
        result_sublabel = min(1, max(0, solution));
        energy_raw_sublabel = value(E);
        valueD = value(D);
        valueS = value(S);
        %{
        valueY = value(Y);
        valueZ_i = value(Z_i);
        valueZ_j = value(Z_j);
        residual_marginalization = 0;
        for edge = 1:num_edges
            Z_1 = valueZ_i(edge, :)';
            Z_2 = valueZ_j(edge, :)';
            Y_edge = squeeze(valueY(edge, :, :));
            
            err_norm = norm(Y_edge - Z_1 * Z_2');
            residual_marginalization = residual_marginalization + err_norm;  
        end
        
        % fprintf('Smoothness term Approximation residual = %f\n', residual_marginalization);
        %}
        
        fprintf('Time: sublabel refinement solve = %f sec\n', time_solve);
        time_elapsed_sublabel = time_first_constraint + time_marginal + ...
            time_smooth_energy + time_solve;
    else
        display('Hmm, something went wrong!');
        sol.info
        yalmiperror(sol.problem)
    end
end

