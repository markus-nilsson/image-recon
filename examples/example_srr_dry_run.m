

% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (1)
    source_fn = '../../data/srr_mdt_data/Serie_03_t1_mprage_tra_ovinkladforneuronavigering.nii.gz';
    I = mdm_nii_read(source_fn);
    
    % Select a slice
    sz = [size(I, 1) size(I, 2)];
    y_hr = reshape(double(I(:,:,60,:)), prod(sz), size(I, 4))';
end

% 2. Downsampling operators
% ----------------------------------------------------------------------
if (1)
    aspect_ratio = 4;
    n_rot = 8;
    theta = linspace(-40, 40, n_rot) / pi * 180;
    
    clear s_list;
    for c = 1:n_rot
        S = imr_op_rotate2d(sz, theta(c), aspect_ratio);
        S_list{c} = op_obj_sample_image(S, sz); %#ok<SAGROW> 
    end
    
end

% 3. Downsample to produce virtual measurement
if (1)
    O = op_obj_sample_images(S_list);
    y_lr = O * y_hr;
end


% 4 Reconstruct and view
if (1)
    clear opts;
    cost.do_mu_update = 0;
    opts.cost = {cost_imfilter(sz, 2, 1e-2)};
    opts.n_iter_admm = 50;

    data = y_lr.add_noise_gaussian(10);

    y_est = srr_est_admm(O, data, opts);

    g = @(x) reshape(x, sz);
    msf_imagesc(cat(1, g(y_est), g(y_hr), 10*abs(g(y_est-y_hr))));
    caxis([0 quantile(y_hr(:), 0.98)]);
end


