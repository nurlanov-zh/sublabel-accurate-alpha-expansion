function [time_elapsed, energy, energy_discretized ...
        ] = Mollenhoff_denoising(L, coef_imresize, lmb, discretize_sublabel)
    %% Truncated quadratic dataterm + TV denoising (see Fig.7) with sublabel lifting
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
    
    fprintf('\n__________\n')
    fprintf('Starting sublabel method by Mollenhoff on GPU...\n')
    % Fixed parameters
    nu = 0.025; % a jump higher than sqrt(0.025) is considered an outlier
    alpha = 25;

    % load image
    im = imread('data/watercastle.jpg');
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

    
    tic;
    gamma = linspace(0, 1, L)';

    if (gamma(2) - gamma(1)) > 2 * sqrt(nu)
        error('ERROR: Not enough labels for dataterm!')
    end

    % compute piecewise quadratic approximation of dataterm
    polya = zeros(L - 1, ny, nx);
    polyb = zeros(L - 1, ny, nx);
    polyc = zeros(L - 1, ny, nx);

    for i=1:(L-1)
        % completely piecewise constant part
        case1 = (gamma(i + 1) <= (im_noisy - sqrt(nu))) | ...
                (gamma(i) >= (im_noisy + sqrt(nu)));

        % completely quadratic part
        case2 = (gamma(i) >= (im_noisy - sqrt(nu))) & ...
                (gamma(i + 1) <= (im_noisy + sqrt(nu)));

        % on left kink
        case3 = (gamma(i) < (im_noisy - sqrt(nu))) & ...
                (gamma(i+1) > (im_noisy - sqrt(nu)));

        % on right kink
        case4 = (gamma(i) < (im_noisy + sqrt(nu))) & ...
                (gamma(i+1) > (im_noisy + sqrt(nu)));

        polya(i, case1) = 0;
        polyb(i, case1) = 0;
        polyc(i, case1) = alpha * nu / 2;

        polya(i, case2) = alpha / 2;
        polyb(i, case2) = -alpha * im_noisy(case2);
        polyc(i, case2) = (alpha / 2) * im_noisy(case2) .^ 2;

        polya(i, case3) = 0;
        polyb(i, case3) = ((alpha / 2) * (gamma(i + 1) - im_noisy(case3)) .^2 - (alpha * nu) / 2) / ...
            (gamma(i + 1) - gamma(i));
        polyc(i, case3) = (alpha * nu / 2) - polyb(i, case3) * gamma(i);

        polya(i, case4) = 0;
        polyb(i, case4) = ((alpha / 2) * (gamma(i) - im_noisy(case4)) .^2 - (alpha * nu) / 2) / ...
            (gamma(i) - gamma(i+1));
        polyc(i, case4) = (alpha * nu / 2) - polyb(i, case4) * gamma(i + 1);
    end


    %% solve problem
    %u_unlifted = sublabel_lifting_quad(polya, polyb, polyc, gamma, lmb);
    
    u_unlifted = sublabel_lifting_quad_l1TV(polya, polyb, polyc, gamma, lmb);
    
    time_elapsed = toc;
    fprintf('Time elapsed Mollenhoff = %f sec\n', time_elapsed);

    %% compute truncated ROF energy
    Kmat = spmat_gradient2d(nx,ny,1);
    [m, n] = size(Kmat);
    grad = Kmat * u_unlifted(:);
    
    % gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
    gradnorms = sum(abs(grad(:)));
    
    energy = 0.5 * alpha * sum(min((u_unlifted(:)-im_noisy(:)).^2, nu)) ...
              + lmb * sum(gradnorms);
    fprintf('Energy Mollenhoff = %f\n', energy);
    
    im_Moll = reshape(u_unlifted(:), [ny, nx]);
    imwrite(im_Moll, ['results/images/sublabel_Mollenhoff_' num2str(L) '.png'])
    imshow(im_Moll);
    
    if discretize_sublabel
        [dis_baseline, ~] = roundtowardvec(u_unlifted, gamma);
        [energy_discretized, compD_dissub, compS_dissub ...
            ] = compute_energy_truncrof(dis_baseline, im_noisy, lmb);
        
        fprintf('Energy discretized Mollenhoff = %f, D = %f, S = %f\n',  ...
            energy_discretized, compD_dissub, compS_dissub);
        im_dis_sublabel = reshape(dis_baseline, [ny, nx]);
        imwrite(im_dis_sublabel, ['results/images/sublabel_Mollenhoff_discretized_' num2str(L) '.png'])
    else
        energy_discretized = -1;
    end
end
