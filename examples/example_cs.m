

% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (1)
    source_fn = '../../data/srr_mdt_data/Serie_03_t1_mprage_tra_ovinkladforneuronavigering.nii.gz';
    data_hr = do_w_from_nii(source_fn);
    data_hr.trim([], [], 90, []);
end

% 2. Compile the forward process, random sampling
if (1)

    % Fourier transform operator
    O_FT = op_obj_image_ft(data_hr.h, data_hr.h, 1);

    % Random sample operator
    rng(1);
    W = rand(size(data_hr.w)) > 0.4; % 40% sampling  
    W(1) = 1; % k-space center frequency

    O_RS = op_obj_image_reweigh(W, W, data_hr.h, data_hr.h, 1);
    
    % Compile the forward process and compute sampled data
    O = op_obj_append(O_FT, O_RS);    

    data = O * data_hr;

end


% 4 Reconstruct and view
if (1)
    clear opts;
    opts.cost = {...
        cost_positivity(0.5), ...
        cost_imfilter_3d(5, 0.04), ...        
        cost_lasso(0.05, 1, 10)};
    opts.n_iter_admm = 20;
    opts.n_iter_cg = 8;

    y_est = srr_est_admm(O, data, opts);

    msf_imagesc(abs(y_est.imreshape()));

end


