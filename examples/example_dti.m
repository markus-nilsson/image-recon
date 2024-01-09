% Purpose:
% 
% Reconstruct the diffusion tensor from a DTI data set, using a 
% log-operator. Warning: Does not work, needs more work.

% 1. Load DTI set
% ----------------------------------------------------------------------
if (1)
    source_fn = fullfile(example_data_path, ...
        '/dti/Serie_06_DTI_30dir_mc.nii.gz');

    data = do_w_from_nii(source_fn);

    data.trim(20:110, 10:115, 30 + (-3:3), []);

    data_log = do_w_image_log(data.w, data.h);

end

% 2. Kernel and log operators
% ----------------------------------------------------------------------
if (1)
    O = op_obj_kernel_dti_cum(data);
    L = op_obj_image_logexp(data.w, data.h, data.h, O.n_k);
    O2 = op_obj_append(O, L);
end

% 3. Reconstruct and view
if (1)
    clear opts;

    opts.cost = {...
        cost_imfilter_3d(1, 0.02)};
    
    opts.n_iter_admm = 10;
    opts.n_iter_cg = 10;

    y_est = srr_est_admm(O2, data_log, opts);

    dt = y_est.new(y_est.w(:,2:7));    
    fa = dt.new( tm_fa(dt.w));
    md = dt.new( tm_md(dt.w));

    data_rec = O2 * y_est;

    msf_clf;
    subplot(2,2,1);
    msf_imagesc(fa.imreshape())
    caxis([0 1]);

    subplot(2,2,2);
    msf_imagesc(md.imreshape())
    caxis([0 3.0]);

    subplot(2,1,2);
    msf_imagesc(cat(1, data.imreshape(), data_rec.imreshape()), 3, [], 1);

end
