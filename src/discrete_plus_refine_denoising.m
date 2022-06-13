function [ time_elapsed_GCO, energy_to_compare_GCO, ...
    time_elapsed_sublabel, energy_to_compare_sublabel, ...
    energy_raw_sublabel, count_non_convex, energy_to_compare_dissublabel ...
    ] = discrete_plus_refine_denoising( L, coef_imresize, lmb, trunc_smooth, ...
    k_smooth, discrete_method, sublabel_method, discretize_sublabel )

    if ~exist('coef_imresize','var') || isempty(coef_imresize)
      coef_imresize=0.5;
    end
    if ~exist('lmb','var') || isempty(lmb)
      lmb=0.6;
    end
    if ~exist('trunc_smooth','var') || isempty(trunc_smooth)
      trunc_smooth=-1;
    end
    if ~exist('k_smooth','var') || isempty(k_smooth)
      k_smooth=1;
    end
    if ~exist('discrete_method','var') || isempty(discrete_method)
      discrete_method='GCO';
    end
    if ~exist('sublabel_method','var') || isempty(sublabel_method)
      sublabel_method='QL';
    end
    if ~exist('discretize_sublabel','var') || isempty(discretize_sublabel)
      discretize_sublabel=true;
    end

    im_noisy = [];
    label_space_double = linspace(0, 1, L)'; % in [0, 1]
    return_neighbors = true;
    
    % 0. Initialize problem
    [data_cost, smooth_cost, neighbor_cost, ~,...
            nx, ny, neighbor_rows, neighbor_cols, im_noisy] = truncrof_init(label_space_double, ...
            return_neighbors, coef_imresize, im_noisy, lmb, trunc_smooth, k_smooth);
        
    fprintf('Image size = %dx%d = %d\n', nx, ny, nx*ny);
    fprintf('Labels num = %d\n', L);

    init_GCO = [];

    %% 1. Run discrete solver
    if strcmp(discrete_method, 'exact')
        % Run the exact solver
        [label_inds_GCO, energy_raw_GCO, time_elapsed_GCO, D, S] = run_multilabel_exact(...
            data_cost, smooth_cost, neighbor_cost, init_GCO);
    elseif strcmp(discrete_method, 'GCO')
        [label_inds_GCO, energy_raw_GCO, time_elapsed_GCO, D, S] = run_GCO(...
            data_cost, smooth_cost, neighbor_cost, init_GCO);
    else
        msg = 'One of the 2 options for `discrete_method` should be chosen: "exact", "GCO".';
        error(msg);
    end
    
    % obtain the values of label indexes
    result_GCO_double = zeros(length(label_inds_GCO), 1);
    for i = 1:length(label_inds_GCO)
        result_GCO_double(i) = label_space_double(label_inds_GCO(i));
    end

    fprintf('Raw energy GCO = %f, D = %f, S = %f\n', energy_raw_GCO, D, S);
    
    [energy_to_compare_GCO, compD_GCO, compS_GCO] = compute_energy_truncrof(...
        result_GCO_double, im_noisy, lmb, trunc_smooth, k_smooth);
    
    fprintf('True energy GCO = %f, D = %f, S = %f\n', energy_to_compare_GCO,...
        compD_GCO, compS_GCO);

    im_GCO = reshape(result_GCO_double, [ny, nx]);
    imwrite(im_GCO, ['results/images/discrete_' discrete_method '_' num2str(L) '.png'])
    %imshow(im_GCO);

    
    %% 2. Run sublabel refinement   
    if strcmp(sublabel_method,'LM') 
        [ result_sublabel_double, time_elapsed_sublabel, energy_raw_sublabel, ...
            subD, subS, count_non_convex ] = sublabel_local_LM(label_inds_GCO, ...
            label_space_double, data_cost, smooth_cost, neighbor_rows, ...
            neighbor_cols);
    elseif strcmp(sublabel_method,'QM')
        [ result_sublabel_double, time_elapsed_sublabel, energy_raw_sublabel, ...
            subD, subS, count_non_convex ] = sublabel_local_QM(label_inds_GCO, ...
            label_space_double, data_cost, smooth_cost, neighbor_rows, ...
            neighbor_cols);
    elseif strcmp(sublabel_method,'QL')
        [ result_sublabel_double, time_elapsed_sublabel, energy_raw_sublabel, ...
            subD, subS, count_non_convex] = sublabel_local_QL(label_inds_GCO, ...
            label_space_double, data_cost, smooth_cost, neighbor_rows, ...
            neighbor_cols);
    else
        msg = 'One of the 3 options for sublabel_method should be chosen: "LM", "QM", "QL".';
        error(msg);
    end

    % Compute energy to compare
    fprintf('Raw energy sublabel = %f, D = %f, S = %f\n', energy_raw_sublabel, subD, subS);
    
    [energy_to_compare_sublabel, compD_sub, compS_sub...
        ] = compute_energy_truncrof(result_sublabel_double, im_noisy, lmb, trunc_smooth, k_smooth);
    
    fprintf('True energy sublabel = %f, D = %f, S = %f\n', energy_to_compare_sublabel, compD_sub, compS_sub);

    im_sublabel = reshape(result_sublabel_double, [ny, nx]);
    imwrite(im_sublabel, ['results/images/sublabel_ours_' sublabel_method ...
        '_' discrete_method '_' num2str(L) '.png'])
    % imshowpair(im_GCO, im_sublabel, 'montage');
    
    %% Discretize sublabel accurate result
    if discretize_sublabel
        [discretized_sublabel, ~] = roundtowardvec(result_sublabel_double, label_space_double);
        [energy_to_compare_dissublabel, compD_dissub, compS_dissub ...
            ] = compute_energy_truncrof(discretized_sublabel, im_noisy, lmb, trunc_smooth, k_smooth);
        fprintf('True energy discretized sublabel = %f, D = %f, S = %f\n',  ...
            energy_to_compare_dissublabel, compD_dissub, compS_dissub);
        im_dis_sublabel = reshape(discretized_sublabel, [ny, nx]);
        imwrite(im_dis_sublabel, ['results/images/discretized_ours_' ...
            sublabel_method '_' discrete_method '_' num2str(L) '.png'])
    else
        energy_to_compare_dissublabel = -1;
    end
end

