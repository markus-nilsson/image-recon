% Purpose: Illustrate SRR using the ADMM framework


% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (1)
    
    source_fn = fullfile(example_data_path(), ...
        'srr_mdt_data/Serie_03_t1_mprage_tra_ovinkladforneuronavigering.nii.gz');

    data_hr = do_w_from_nii(source_fn);
    data_hr.trim([], [], 60, []); % select one slice only for now


end

% 2. Downsampling operators
%    (Not a perfect representation of real data, but it works for now)
% ----------------------------------------------------------------------
if (1)
    aspect_ratio = 4;
    n_rot = 8;
    theta = linspace(-40, 40, n_rot) / 180 * pi;

    h0 = data_hr.h;
    h0.srow_x(1:3) = [1 0 0];
    h0.srow_y(1:3) = [0 1 0];
    h0.srow_z(1:3) = [0 0 1];
    h0.qform_code = 0;   
    
    clear s_list;
    for c = 1:n_rot

        t = theta(c);
        h = h0;

        r = [h.srow_x(4) h.srow_y(4) h.srow_z(4)];

        M = [...
            cos(t) sin(t) 0;
            -sin(t) cos(t) 0;
            0 0 1];

        M = M * diag([aspect_ratio 1 1 ]);

        r0 = M * r';

        h.srow_x(1:4) = [M(1,:)  r0(1)];
        h.srow_y(1:4) = [M(2,:)  r0(2)];
        h.srow_z(1:4) = [M(3,:)  r0(3)];
  
        S = op_obj_image_h2l(h, h0, 1);
        
        S_list{c} = S;
    end
    
end

% 3. Downsample to produce virtual measurement
if (1)
    O = op_obj_image_stack(S_list);
    data = O * data_hr;
end


% 4 Reconstruct and view
if (1)
    clear opts;
    cost.do_mu_update = 0;
    opts.cost = {cost_imfilter_3d(1, 0.05)};

    opts.n_iter_admm = 15;

    y_est = srr_est_admm(O, data, opts);

    msf_imagesc((cat(1, data_hr.imreshape(), y_est.imreshape()))); colorbar
end


