function [data_corr, p] = mec_run(data, opt)

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

p.global = mec_params_global_orders(data, opt.ec_order);
p.volume = mec_params_volume_rigid(data);
p.slice  = mec_params_slice(data);

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
        p.global = p.global.optimize(data, data_ref, points);
        p.global.plot();
    end

    % Volume fitting
    p.volume.prior_param = p.global;
    
    if (~opt.volume.skip)
        p.volume = p.volume.optimize(data, data_ref, points);
        p.volume.plot();
    end

    % Slice fitting
    p.slice.prior_param = p.volume;
    
    if (~opt.slice.skip)
        p.slice = p.slice.optimize(data, data_ref, points);
        p.slice.plot();
    end


    % Apply to all data
    data_corr = mec_points_all(data).apply(mec_data(data, 'cubic'), p.slice);

    % Plot it
    if (0)

        figure(my_fig('mec_run'));
        msf_clf;
        ind = data.xps.b > (0.9 * max(data.xps.b));
        A = std(data.trim_copy([],[],[],ind).imreshape(), [], 4);
        B = std(data_corr.trim_copy([],[],[],ind).imreshape(), [], 4);
        C = std(data_ref.trim_copy([],[],[],ind).imreshape(), [], 4);

        C = C * (C(:) \ A(:));

        k_list = round([0.15 0.5 0.85] * data.dim(3));
        J = [];
        for k = k_list
            J = cat(2, J, cat(1, A(:,:,k), B(:,:,k), C(:,:,k)));
        end
        msf_imagesc(J);

       
        figure(my_fig('mec_global'));
        msf_clf;
        p.plot_global();

        figure(my_fig('mec_volume'));
        msf_clf;
        p.plot_volume();

        % figure(4); msf_clf;
        % p.plot_slice();
        % 
        % 
        % figure(f);

        pause(0.05);
    end

end

