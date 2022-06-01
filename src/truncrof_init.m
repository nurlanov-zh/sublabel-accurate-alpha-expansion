function [data_cost, smooth_cost, neighbor_cost, time_elapsed,...
    nx, ny, index_rows, index_cols, im_noisy] = truncrof_init(label_space, ...
    return_neighbors, coef_imresize, im_noisy, lmb, trunc_smooth, k_smooth)
    
    % denoising problem with truncated ROF model
    if ~exist('return_neighbors','var') || isempty(return_neighbors)
      return_neighbors=true;
    end
    if ~exist('coef_imresize','var') || isempty(coef_imresize)
      coef_imresize=0.5;
    end
    if ~exist('im_noisy','var') || isempty(im_noisy)
      im_noisy=[];
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
    
    % Fixed parameters from previous work by Mollenhoff (2016).
    nu = 0.025;
    alpha = 25.0;
    
    % fprintf('\nStarting truncated RoF problem construction...\n');
    %% Truncated quadratic dataterm + TV denoising
    if isempty(im_noisy)
        rng(42);
        % load image
        im_path = 'watercastle.jpg';
        im = imread(im_path);
        im = double(imresize(im, coef_imresize)) / 255;
        [ny, nx] = size(im);

        N = ny * nx;

        % noise parameters
        noise_sigma = 0.05; % standard deviation of gaussian noise
        noise_sp = 0.25;    % percentage of salt&pepper noise

        % add gaussian noise
        im_noisy = im + noise_sigma * randn(ny, nx, 1);

        % add salt and pepper noise
        im_noisy = im_noisy(:);
        perm = randperm(N);
        num_sp = round(N *noise_sp * 0.5);
        im_noisy(perm(1:num_sp)) = 1;
        im_noisy(perm(num_sp+1:2*num_sp)) = 0;
        im_noisy = min(max(im_noisy, 0), 1);
        im_noisy = reshape(im_noisy, [ny, nx]);
    end
    
    [L, ~] = size(label_space);
    [ny, nx] = size(im_noisy);
    N = ny * nx;
    % fprintf('Image size = %dx%d = %d\n', nx, ny, nx*ny);
    % fprintf('Labels num = %d\n', L);
    
    tic;
    %% Set data cost
    data_cost = zeros(N, L);
    
    for i=1:L
        data_cost(:, i) = 0.5 * alpha ...
                          * min((label_space(i) - im_noisy(:)).^2, nu);
    end
    
    data_cost = data_cost';                       % Should be (NumLabels, NumSites)
    
    %% Set smoothness cost (L1 Total Variation) 
    smooth_cost = zeros(L, L);
    if trunc_smooth > 0
        % truncated smoothness cost
        for i=1:L
            for j = 1:L
                smooth_cost(i, j) = lmb * min(trunc_smooth, (abs(label_space(i) - label_space(j))) .^ k_smooth);
            end
        end
    else
        % L1-TV cost
        for i=1:L
            for j = 1:L
                smooth_cost(i, j) = lmb * abs(label_space(i) - label_space(j));
            end
        end
    end
    
    %% Set neighborhood
    if return_neighbors
        index_rows = zeros(1, 2 * N - nx - ny);
        index_cols = zeros(1, 2 * N - nx - ny);
        values = zeros(1, 2 * N - nx - ny);

        idx = 1;
        for i=1:N
            % Add lower neighbor
            if (mod(i, ny) ~= 0 && i + 1 <= N)
                index_rows(idx) = i;
                index_cols(idx) = i + 1;
                values(idx) = 1.0;
                idx = idx + 1;
            end
            % Add neighbor to the right
            if (i + ny <= N)
                index_rows(idx) = i;
                index_cols(idx) = i + ny;
                values(idx) = 1.0;
                idx = idx + 1;
            end
        end
        neighbor_cost = sparse(index_rows, index_cols, values, N, N);
    else
        neighbor_cost = [];
        index_rows = [];
        index_cols = [];
    end
    
    time_elapsed = toc;
    fprintf('\n________________________________________\n');
    fprintf('Time: truncated ROF construction = %f sec\n', time_elapsed);
    
    %{
    [energy, compD, compS ...
        ] = compute_energy_truncrof(im, im_noisy, lmb);
    fprintf('Energy orig = %f, D = %f, S = %f \n', energy, compD, compS);
    %}
end