function [result_exact, energy_exact, time_elapsed, D, S] = run_multilabel_exact(data_cost,...
    smooth_cost, neighbor_cost, init)

    % data_cost: size (L, N)
    % smooth_cost: size (L, L);
    % neighbor_cost: size (N, N);
    % init is initial labelling, currently not used

    % we always work here with regular TV-L1 because of submodularity
    % constraint
    
    if ~(exist('MultiLabelSubModular_mex.mexa64', 'file') == 3)
        install_mlsm;
    end
        

    tic;
    data_cost = data_cost';
    nl = size(data_cost,2);
    
    % smooth_cost = smooth_cost * (nl - 1);
    data_cost = data_cost * (nl - 1);
    
    % scale_factor = 1e6;
    % data_cost_int = int32(round(scale_factor * data_cost));
    
    lmb = max(max(smooth_cost));
    W = spfun(@(x) lmb, neighbor_cost);  
    % W_int = spfun(@(x) round(scale_factor * lmb), neighbor_cost);
    
    Vi = spfun(@(x) 1, neighbor_cost);
    Vm = abs( bsxfun(@minus, 1:nl, (1:nl)') );
    % Vm_int = int32(round(Vm));
    
    [x, e] = MultiLabelSubModular(data_cost, W, Vi, Vm);
    time_elapsed = toc;
    result_exact = x';          % transpose to (N, 1)
    e = e / (nl - 1);
    energy_exact = e(1);
    D = e(2);
    S = e(3);
end

