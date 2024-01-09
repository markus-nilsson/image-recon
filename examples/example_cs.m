% Purpose: Compressed sensing example


% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (1)

    bp = example_data_path();

    source_fn = fullfile(bp, 'srr_mdt_data/Serie_03_t1_mprage_tra_ovinkladforneuronavigering.nii.gz');
    data_hr = do_w_from_nii(source_fn);
    data_hr.trim([], [], 90, []); % select one slice only
end

% 2. Compile the forward process, random sampling
if (1)
    
    % Deterministic randomness
    rng(10);

    % Fourier transform operator
    O_FT = op_obj_image_ft(data_hr.h, data_hr.h, 1);

    % Construct random sampling operator (i.e. set k-space points to zero)
    rng(1);
    W = rand(size(data_hr.w)) > 0.4; % 60% sampling  
    W(1) = 1; % k-space center frequency

    O_RS = op_obj_image_reweigh(W, W, data_hr.h, data_hr.h, 1);
    
    % Compile the forward process and compute sampled data
    O = op_obj_append(O_FT, O_RS);    

    data = O * data_hr;

end


% 4 Reconstruct and view
if (1)

    % Try changing these arguments – view the functions to know which is
    % which. The value of the arguments affect the performance. 
    clear opts;
    opts.cost = {...
        cost_positivity(0.1, 1), ...
        cost_imfilter_3d(1, 0.01), ...  % 1 - median filtering      
        cost_lasso(2, 1, 10)};
    
    opts.n_iter_admm = 50;
    opts.n_iter_cg = 10;

    y_est = srr_est_admm(O, data, opts);

    msf_imagesc(real(y_est.imreshape()));

end


