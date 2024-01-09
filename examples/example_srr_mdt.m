

% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (1)
    % Build data object
    bp = '../../data/srr_mdt_data/';
    fns = msf_find_fns(bp, '*pa.nii.gz');
    data_lr = data_obj_from_nii_fns(fns);
    data_lr.trim([], 70, [], 3);    
end

% 2. Define sampling operators
if (1)

    % Create an empty hr object for later use, from first lr object
    data_hr = data_obj_image_volume.new_hr_from_lr(data_lr.data_obj{1});

    O = op_obj_image_samplers(data_lr.h, data_hr.h);

    % Show naive upsamplers
    g = @(x) x.imreshape();
    for k = 1:10
        for c = 1:6
            msf_imagesc(g(O.apply_adjoint(data_lr, c)),2), 
            pause(0.05); 
        end
    end
end

% 3. Reconstruct and view
if (0)
    clear opts;
    opts.cost = {cost_imfilter(2, 0.02)};
    opts.n_iter_admm = 15;

    y_est = srr_est_admm(O, data_lr, opts);

    msf_imagesc(y_est.imreshape(),2);
    caxis([0 20]);
    colormap gray;
end


% 4. Find the right level for the regularizatino
if (0)

    p = linspace(-4, 2, 100);

    rd = zeros(size(p)) + NaN;
    rc = zeros(size(p)) + NaN;
    for c = 1:numel(p)

        clear opts;
        opts.present = 1;
        cost = cost_imfilter(sz, 2, 10^p(c));
        cost.do_mu_update = 0;
        opts.cost = {cost};
        opts.n_iter_admm = 50;

        y_est = srr_est_admm(O, data, opts);

        rd(c) = norm( O * y_est - data );

        [f,b] = cost.do_iter(y_est);
        rc(c) = norm(f(y_est) - b);

        g = @(x) reshape(x, sz);
        

        subplot(2,2,1);
        msf_imagesc(cat(1, g(y_est)));

        subplot(2,2,2);
        msf_imagesc(g(cost.u));

        subplot(2,2,3);
        semilogy(p, rd);

        subplot(2,2,4);
        semilogy(p, rc);

        pause(0.05);
    end


end

