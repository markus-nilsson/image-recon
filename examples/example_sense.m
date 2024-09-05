% Purpose: Test the approach for parallell imaging – it works for a
% paralell imaging factor of 2, but not higher, at present

% 1. Load high-resolution images
% ----------------------------------------------------------------------
if (1)
    source_fn = fullfile(example_data_path, ...
        'srr_mdt_data/Serie_03_t1_mprage_tra_ovinkladforneuronavigering.nii.gz');
    data_hr = do_w_from_nii(source_fn);
    data_hr.trim([], [], 90, []); % select one slice only for now

end


% 2. Make our own coil sensitivity maps
if (1)

    n = 32;
    x = (1:data_hr.h.dim(2))';
    y = (1:data_hr.h.dim(3));

    r = 1 * (mean(x) + mean(y)) / 2;

    x0 = -2;
    y0 = 20;

    theta = linspace(0, 2*pi, n+1);   

    % mask from e.g. lr data
    M = mio_mask_expand(data_hr.imreshape() > 30);
    M = M(:,:,:,1);
%     M(M(:) >= 0) = 1;

    O = cell(1, n); L = cell(1, n);
    for c = 1:n

        r0 = [cos(theta(c)) * r + mean(x) - x0 sin(theta(c)) * r + mean(y) - y0];

        n = r0 - [mean(x) mean(y)]; % normal
        n = n / norm(n);

        coord = cat(3, ...
            (x  - r0(1)) * ones(size(x))', ...
            ((y' - r0(2)) * ones(size(y'))')');

        coord = coord ./ repmat( sqrt(sum(coord.^2,3)), [ 1 1 2] );

        area = abs((coord(:,:,1) * n(1) + coord(:,:,2) * n(2)));

%         area = area;

        msf_imagesc(area);
        caxis([0 1]);
        pause(0.05);

        w = (x - r0(1)).^2 + (y - r0(2)).^2;
        w = w / r^2;
        w = w * 20;
        w(w < 1) = 1;
        p = 1.0;
        w = 1./w.^(p);

        w = w .* area;

        W = reshape(w, [numel(w) 1]);
        WT = W;

        W(M(:) == 0) = 0;
        WT(M(:) == 0) = 0;



        O{c} = op_obj_image_reweigh(W, WT, data_hr.h, data_hr.h, 1);
        L{c} = [ 1 ]; 

    end
    
    O_COIL =  op_obj_image_stack(O, L);
end

% Compile the forward process
if (1)

    % Fourier transform operator
    O_FT = op_obj_image_ft(data_hr.h, data_hr.h, 2);

    % Undersampling operator
    PI_factor = 5;
    w = (x > 0) & (mod(y, PI_factor) == 1);
    W = reshape(w, [numel(w) 1]);
    WT = W;   

    O_US = op_obj_image_reweigh(W, WT, data_hr.h, data_hr.h, 1);
    
    % Compile the forward process and compute sampled data
    O = op_obj_append(O_FT, O_US);    

    O2 = op_obj_append(O_COIL, O);
    
    data = O2 * data_hr;

end

% Run sum-of-squares recon
if (0)

    for c = 1:n
        tmp = O_COIL.O_list{c} * (O_FT' * data.data_obj{c});

        if (c == 1), tmp2 = zeros(size(tmp.w)); tmp3 = tmp2; end
        tmp2 = tmp2 + tmp.w.^2;% .* O_COIL.O_list{c}.W;
        tmp3 = tmp3 + O_COIL.O_list{c}.W;
    end
    tmp2 = tmp.new( (abs(sqrt(tmp2)) ./ (tmp3.^2 + eps)));

    tmp2 = tmp2 * (tmp2.w(:) \ data_hr.w(:));

    msf_imagesc(tmp2.imreshape());
    

    error('stop');
end


% 4 Reconstruct and view
%
% https://github.com/marcsous/parallel/blob/8b6f18107eeba424eb662ff5d4f52c82e2407eb2/homodyne.m#L171
% check code above
% adds some extra regularizers, e.g., a L2 regularizer on the total
% solution
% 
% and 
%
% something penalizing the imaginary term overall (e.g. adding a term
% 1i * imag(x) 
%
% this does not have to be an ADMM-like term
if (1)
    clear opts;  

    opts.n_iter_admm = 100; % admm not needed here
    opts.n_iter_cg = 25;

    opts.cg_display_ind = 1;

    opts.cost = {...            
        cost_imfilter_3d(6, 0.001), ...
        };    


%         cost_positivity(0.01), ...
%         cost_sense(0.0), ...
    
%     opts.init_x = @(x) tmp2;

    [y_est,r,cgr] = srr_est_admm(O2, data, opts);


    srr_plot_convergence(cgr);

    msf_clf;
    msf_imagesc(abs(y_est.imreshape()));
    colorbar;

end


