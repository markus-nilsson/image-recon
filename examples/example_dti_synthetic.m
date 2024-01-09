% Purpose: Reconstruct DTI data from multiple independent rotations
% 
error('Not compatible with this verision of the framework, yet'); 

% 1. Load DTI set
% ----------------------------------------------------------------------
if (1)
    bp = fullfile(example_data_path, 'dti');
    source_fn = fullfile(bp, 'Serie_06_DTI_30dir_mc.nii.gz');
    I = mdm_nii_read(source_fn);

    xps = mdm_xps_from_gdir(fullfile(bp, 'Serie_06_DTI_30dir_mc_gdir.txt'));
    
    % Select a slice
    sz = [size(I, 1) size(I, 2)];
    y_hr = reshape(double(I(:,:,30,:)), prod(sz), size(I, 4))';
end

% 2. Downsampling operators
% ----------------------------------------------------------------------
if (1)
    aspect_ratio = 3;
    n_rot = 10;
    theta = linspace(-90, 90, n_rot) / pi * 180;
    
    S_list = {};
    for c = 1:n_rot
        
        R = op_obj_image_h2l.tmp_rotate2d(sz, theta(c), aspect_ratio);

        k_list = [1 1 + (1:3) + (c-1)*3];
        for k = 1:numel(k_list)
            S = op_obj(R, sz);
            S.k = k_list(k);
            S_list{end+1} = S;
        end
    end
    
end

% 3. Downsample to produce virtual measurement
if (1)
    O = op_obj_image_stack(S_list);
    y_lr = O * y_hr;
end

% 4. Add a DTI kernel 
if (1)
    K = op_obj_kernel_dti(xps);
    J = op_obj_append(K, O);
end


% 5. Reconstruct and view
if (1)
    clear opts;
    opts.cost = {...
        cost_imfilter(sz, 2, 1e1, 0), ...
        cost_dict(K, 1e-4)};
    opts.n_iter_admm = 25;

    y_est = srr_est_admm(J, y_lr, opts);

    y_rec = K * y_est;

    msf_clf;
    for c = 1:size(y_hr,1)
        g = @(x) reshape(x, sz);
        msf_imagesc(cat(1, g(y_rec(c,:)), g(y_hr(c,:)), 7*abs(g(y_rec(c,:)-y_hr(c,:)))));
        caxis([0 quantile(y_hr(c,:), 0.99)]);
        pause(0.15);
    end
end

% 6. Grab the diffusion metrics (works in progress)
if (0)
    tmp = K.x_to_m(y_est);
end 

