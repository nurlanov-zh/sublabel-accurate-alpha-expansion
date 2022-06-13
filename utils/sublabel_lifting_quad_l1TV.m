function u_unlifted = sublabel_lifting_quad_l1TV(polya, polyb, polyc, gamma, lmb)
% this is an implementation of the proposed sublabel lifting method
% for a dataterm that is given by a piecewise quadratic energy.
    
    [k, ny, nx] = size(polya);
    L = size(gamma, 1);
    N = ny*nx;

    if k ~= (L-1)
        fprintf(['Label space and size of polynomials doesnt fit ' ...
                 'together (%d,%d)\n'], k, L-1);
        return;
    end
    
    polya = polya(:);
    polyb = polyb(:);
    polyc = polyc(:);

    gamma_vec_start = repmat(gamma(1:end-1), [N, 1]);
    gamma_vec_end = repmat(gamma(2:end), [N, 1]);

    lmb_scaled = lmb * (gamma(2) - gamma(1));

    %% setup problem

    % primal variables
    u = prost.variable(N*k);
    s = prost.variable(N*k);
    
    % dual variables
    vz = prost.variable(2*N*k);
    q = prost.variable(N);
    p = prost.variable(2*N*k);
    
    v = prost.sub_variable(vz, N*k);
    z = prost.sub_variable(vz, N*k);

    problem = prost.min_max_problem( {u,s}, {vz, q, p} );

    problem.add_function(vz, prost.function.transform( @(idx, count) ...
                                                       prox_sum_ind_epi_conjquad_1d(idx, count / 2, ...
                                                          false, polya, ...
                                                          polyb, polyc, ...
                                                          gamma_vec_start, ...
                                                          gamma_vec_end), ...
                                                       1 / (gamma(2) - gamma(1)), 0, 1, 0, 0));
    
    problem.add_function(q, prost.function.sum_1d('zero', 1, 0, 1, 1, 0));
    
    %problem.add_function(p, prost.function.sum_norm2(2, false, 'ind_leq0', 1/lmb_scaled, 1, ...
    %                                                 1, 0, 0));
    
    
    
    problem.add_function(p, prost.function.sum_1d('ind_box01', 0.5 / lmb_scaled, ...
        -0.5, 1, 0, 0));
    
    %problem.add_function(p, prost.function.sum_norm2(1, false, 'ind_leq0', 1/lmb_scaled, 1, ...
     %   1, 0, 0));
    
    
    

    problem.add_dual_pair(u, v, prost.block.identity());
    problem.add_dual_pair(s, v, @(row, col, nrows, ncols) ...
                            block_dataterm_sublabel(row, col, nx, ny, ...
                                                    L, gamma(1), ...
                                                    gamma(end)));
    problem.add_dual_pair(s, z, prost.block.identity());
    problem.add_dual_pair(s, q, prost.block.sparse(kron(speye(N), (gamma(1:end-1) - ...
                                                      gamma(2:end))')));
    problem.add_dual_pair(u, p, prost.block.gradient2d(nx, ny, k, true));
    

    %% create backend
    backend = prost.backend.pdhg(...
        'tau0', 10, ...
        'sigma0', 0.1, ...
        'stepsize', 'boyd');

    %% specify solver options
    opts = prost.options(...
        'max_iters', 25000, ...
        'num_cback_calls', 25, ...
        'tol_rel_primal', 1e-6, ...
        'tol_abs_primal', 1e-6, ...
        'tol_rel_dual', 1e-6, ...
        'tol_abs_dual', 1e-6);
    
    %% solve problem
    solution = prost.solve(problem, backend, opts);

    %% obtain result via layer-cake
    p_val = p.val;
    u_volume = reshape(u.val, [k, N]);
    u_unlifted = sum(u_volume, 1) * (gamma(2) - gamma(1));
end

