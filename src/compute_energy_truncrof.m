function [ energy_rof, data_term, smoothness_term ...
    ] = compute_energy_truncrof(result, im_noisy, lmb, trunc_smooth, k_smooth)
    
    if ~exist('trunc_smooth','var') || isempty(trunc_smooth)
      trunc_smooth=-1;
    end
    if ~exist('k_smooth','var') || isempty(k_smooth)
      k_smooth=1;
    end
    
    % Fixed parameters from prefivous work by Mollenhoff (2016)
    nu = 0.025;
    alpha = 25.0;

    [ny, nx] = size(im_noisy);
    Kmat = spmat_gradient2d(nx,ny,1);
    [~, n] = size(Kmat);
    grad = Kmat * result(:);
    
    % gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
    if trunc_smooth > 0
        gradnorms = sum(min(trunc_smooth, (abs(grad(:))).^k_smooth));
    else
        gradnorms = sum(abs(grad(:)));
    end
    
    data_term = 0.5 * alpha * sum(min((result(:)-im_noisy(:)).^2, nu));
    smoothness_term = lmb * sum(gradnorms);
    energy_rof = data_term + smoothness_term;
end

