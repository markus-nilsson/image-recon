

% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (0)
    source_fn = '../../data/srr_mdt_data/Serie_03_t1_mprage_tra_ovinkladforneuronavigering.nii.gz';
    I = mdm_nii_read(source_fn);
    
    % Select a slice
    I = double(I(1:2:end,1:2:end,60));
    sz = size(I);
    y_hr = reshape(I, prod(sz), 1)';
end

% 2. K-space sampling operator
% ----------------------------------------------------------------------
if (0)

    x = (1:sz(1)); x = x - mean(x);
    y = (1:sz(2))';y = y - mean(y);

    N = numel(x);

    S = zeros(numel(I), numel(I));
    for i = 1:sz(1)
        for j = 1:sz(2)
            F = exp(-2i*pi * (x(i)*x + y(j)*y) / N);
            S(:, (i-1)*sz(1) + j) = F(:);
        end
    end

    % Eliminate some 20% of the samples randomly
    S2 = S(randperm(size(S,1), round(0.8*size(S,1))), :);

    P = op_obj_sample_image(S2, sz);
    O = op_obj_sample_images({P});
    
end

% 3. Sample to produce virtual measurement
if (1)
    y_lr = O * y_hr;
end

% 4 Reconstruct and view -- does not work -- need sparsifying space here!
if (1)
    clear opts;
    opts.present = 1;
    cost = cost_imfilter(sz, 4, 1e-2);
    cost.do_mu_update = 0;
    opts.cost = {cost};
    opts.n_iter_admm = 15;

    data = y_lr;

    y_est = srr_est_admm(O, data, opts);


    g = @(x) reshape(x, sz);
    msf_imagesc(cat(1, g(y_est), g(y_hr), 10*abs(g(y_est-y_hr))));
    caxis([0 quantile(y_hr(:), 0.98)]);
end


