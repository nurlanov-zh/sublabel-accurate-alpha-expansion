function [result_gco, energy_gco, time_elapsed, D, S] = run_GCO(data_cost,...
    smooth_cost, neighbor_cost, init)

    % GCO_BuildLib(struct('Debug',1,'EnergyType','double'));  % too slow!    
    % GCO_BuildLib();
    
    fprintf('\nStarting GCO...\n');
    GCO_LoadLib;  
    
    [L, N] = size(data_cost);
    tic;

    % Create GCO problem
    h = GCO_Create(N, L);             % Create new object with NumSites=N, NumLabels=L

    %% Set data cost
    scale_factor = 1e6;
    
    data_cost_int = int32(round(scale_factor * data_cost));
    GCO_SetDataCost(h, data_cost_int);

    %% Set smoothness cost
    
    smooth_cost_int = int32(round(scale_factor * smooth_cost));
    GCO_SetSmoothCost(h, smooth_cost_int);  

    %% Set neighborhood
    GCO_SetNeighbors(h, neighbor_cost);
    
    if ~isempty(init)
        GCO_SetLabeling(h, init);
    end
    
    %% Compute optimal labeling via alpha-expansion 
    GCO_Expansion(h);
    result_gco = GCO_GetLabeling(h);
    time_elapsed = toc;
    fprintf('Time: GCO solve = %f sec\n', time_elapsed);

    [energy_gco, D, S, ~] = GCO_ComputeEnergy(h);
    energy_gco = double(energy_gco) / scale_factor;
    D = double(D) / scale_factor;
    S = double(S) / scale_factor;

    GCO_Delete(h);                   % Delete the GCoptimization object when finished
end

