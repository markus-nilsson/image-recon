function [data_corr, p] = mec_run2(data, opt)

if (nargin < 2), opt.present = 1; end

opt = msf_ensure_field(opt, 'n_iter', 3);
opt = msf_ensure_field(opt, 'ec_order', 2);

% data is an image vector holding an xps too

% Loop over these steps

% 1. Predict reference (could be b0, or 4d volume)

% 2. Run all or any: global, volume (rigid), slice (intra-volume motion)

% 3. Warp all volumes

% Treat these as different optimization problems, but apply them all in
% the final warp

data_corr = data;

data = mec_data(data, 'linear');

p = mec_params_composite(data, opt.ec_order);

for c_iter = 1:opt.n_iter

    % Predict reference
    switch (c_iter)

        case 1
            opt.ref.c_case = 2; % shell-wise reference
            opt.ref.lambda = 0.15;    

            opt.global.skip = 0;
            opt.volume.skip = 1;
            opt.slice.skip = 1;

            % Optimization points
            mask = data.mask_for_mec();
            points = mec_points_some(data, mask, 0.05);
            

        case 2
            opt.ref.c_case = 1; 
            opt.ref.lambda = 0.10;     

            opt.volume.skip = 0;

            % Optimization points
            mask = data.mask_from_r(data_corr, 0.1);
            points = mec_points_some(data, mask, 1);
           

        case 3 % redo with improved reference
            opt.ref.lambda = 0.07;     
            
        case 4 % enable slice
            opt.slice.skip = 0;
    end

    data_ref = mec_data(mec_ref_predict(data_corr, opt.ref));

    % Global fitting
    if (~opt.global.skip)
        p = p.optimize('global', data, data_ref, points);
        p.p_global.plot();
    end

    % Volume fitting  
    if (~opt.volume.skip)
        p = p.optimize('volume', data, data_ref, points);
        p.p_volume.plot();
    end

    % Slice fitting
    if (~opt.slice.skip)
        p = p.optimize('slice', data, data_ref, points);
        p.p_slice.plot();
    end

    % Apply to all data
    data_corr = mec_points_all(data).apply(mec_data(data, 'cubic'), p);

end

