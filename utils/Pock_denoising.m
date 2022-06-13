function [time_elapsed, energy, energy_discretized ...
        ] = Pock_denoising(L, coef_imresize, lmb, discretize_sublabel)
    %% Truncated quadratic dataterm + TV denoising (see Fig.7) with baseline lifting
    rng(42);
    if ~exist('coef_imresize','var') || isempty(coef_imresize)
      coef_imresize=0.5;
    end
    if ~exist('lmb','var') || isempty(lmb)
      lmb=0.6;
    end
    if ~exist('discretize_sublabel','var') || isempty(discretize_sublabel)
      discretize_sublabel=true;
    end
    
    % Fixed parameters
    nu = 0.025; % a jump higher than sqrt(0.025) is considered an outlier
    alpha = 25;

    % load image
    im = imread('data/watercastle.jpg');
    im = double(imresize(im, coef_imresize)) / 255;
    % imwrite(im, 'watercastle_05.jpg');
    [ny, nx] = size(im);
    N = ny * nx;

    fprintf('\n__________\n')
    fprintf('Starting sublabel method by Pock on GPU...\n')
    fprintf('Image size: %d x %d = %d pixels\n', nx, ny, N);
    fprintf('Labels num: %d\n', L);

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
    % imwrite(im_noisy, 'watercastle_noisy_05.jpg');

    tic;
    gamma = linspace(0, 1, L)';

    % create cost volume
    cost_volume = zeros(ny, nx, L);
    for i=1:L
        cost_volume(:, :, i) = (alpha / 2) * min((gamma(i) - im_noisy).^2, nu);
    end

    %% solve problem and display result
    u_unlifted = baseline_lifting_l1TV(cost_volume, gamma, lmb);

    time_elapsed = toc;
    fprintf('Time elapsed Pock = %f sec\n', time_elapsed);

    %% compute truncated ROF energy
    Kmat = spmat_gradient2d(nx,ny,1);
    [m, n] = size(Kmat);
    
    grad = Kmat * u_unlifted(:);
    
    % gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
    gradnorms = sum(abs(grad(:)));
    energy = 0.5 * alpha * sum(min((u_unlifted(:)-im_noisy(:)).^2, nu)) ...
              + lmb * sum(gradnorms);
    fprintf('Energy Pock = %f\n', energy);
    
    
    im_Pock = reshape(u_unlifted, [ny, nx]);
    imwrite(im_Pock, ['results/images/sublabel_Pock_' num2str(L) '.png'])
    
    imshow(im_Pock);
    
    if discretize_sublabel
        [dis_baseline, ~] = roundtowardvec(u_unlifted, gamma);
        [energy_discretized, compD_dissub, compS_dissub ...
            ] = compute_energy_truncrof(dis_baseline, im_noisy, lmb);
        fprintf('Energy discretized sublabel Pock = %f, D = %f, S = %f\n',  ...
            energy_discretized, compD_dissub, compS_dissub);
        im_dis_sublabel = reshape(dis_baseline, [ny, nx]);
        imwrite(im_dis_sublabel, ['results/images/sublabel_Pock_discretized_' num2str(L) '.png'])
    else
        energy_discretized = -1;
    end
end
