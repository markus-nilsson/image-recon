

% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (1)
    source_fn = '../../data/srr_mdt_data/Serie_03_t1_mprage_tra_ovinkladforneuronavigering.nii.gz';
    data_hr = do_w_from_nii(source_fn);
    data_hr.trim([], [], 90, []);
end

% 2. Make coil sensitivity maps
if (1)

    n = 8;
    x = (1:data_hr.h.dim(2))';
    y = (1:data_hr.h.dim(3));

    r = 1.5 * (mean(x) + mean(y)) / 2;

    theta = linspace(0, 2*pi, n+1);   

    O = cell(1, n); L = cell(1, n);
    for c = 1:n

        r0 = [cos(theta(c)) * r + mean(x) sin(theta(c)) * r + mean(y)];
        w = (x - r0(1)).^2 + (y - r0(2)).^2;
        w = w / r^2;
        w = w * 20;
        w(w < 1) = 1;
        w = 1./w;

        W = reshape(w, [numel(w) 1]);
        WT = 1./W;

        O{c} = op_obj_image_reweigh(W, WT, data_hr.h, data_hr.h, 1);
        L{c} = 1;

    end
    
    O_COIL =  op_obj_image_stack(O, L);
end

% Compile the forward process
if (1)

    % Fourier transform operator
    O_FT = op_obj_image_ft(data_hr.h, data_hr.h, 1);

    % Undersampling operator
    w = (x > 0) & (mod(y, 2) == 1);
    W = reshape(w, [numel(w) 1]);
    WT = W;   

    O_US = op_obj_image_reweigh(W, WT, data_hr.h, data_hr.h, 1);
    
    % Compile the forward process and compute sampled data
    O = op_obj_append(O_FT, O_US);    

    O2 = op_obj_append(O_COIL, O);
    
    data = O2 * data_hr;

end


% 4 Reconstruct and view
if (1)
    clear opts;

    opts.n_iter_admm = 5;

    y_est = srr_est_admm(O2, data, opts);

    msf_imagesc(abs(y_est.imreshape()));

end


