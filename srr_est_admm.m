function [x,r,cgr] = srr_est_admm(O, data, opts)
% function x = smr_mrf_admm_recon(O, data, opts)

% Manage  input...
if (nargin < 3), opts = []; end
opts = default_opts(opts, O);

% ADMM Iterations
r = zeros(1, opts.n_iter_admm);
for j = 1:opts.n_iter_admm
    disp(j)

    % build the problem; begin with data consistency term
    % min(x) ||A*x - y||_2^2, here y = data
    % solution when A^T * A * x = A^T * y

    if (j == 1) % Init: y = A * x --> data = E * x (n,1 = n,m x m,1)
        bp = O' * data;
        x = opts.init_x();
    end

    b = bp; % A' * y precomputed
    f = @(x) O' * (O * x);

    % solve secondary problems & update the problem
    for c = 1:numel(opts.cost)
        [fn,bn] = opts.cost{c}.do_iter(x);
        f = @(x) f(x) + fn(x);
        b = b + bn;
    end

    % solve main problem; update x and compute the residual
    x_old = x;
    [x,cgr{j}] = srr_conj_grad(f, b, x, opts);
    r(j) = norm(O * x - data);

    if (j > 2) && (abs(r(j)) > min(abs(r(1:(j-1)))))
        % x = x_old;
        % if (j > 20), break; end
    end

    if (j > 20)
        1;
    end


    subplot(2,2,1);
    srr_plot_convergence(cgr);

    subplot(2,2,3);
    plot(r);
    
    subplot(1,2,2);
    msf_imagesc(x.imreshape());
    pause(0.05);

end

end

function opts = default_opts(opts, O)

    function opts = f(opts, fn, value)
        if (~isfield(opts, fn)), opts.(fn) = value; end
    end

opts = f(opts, 'init_x', @() O.init_x());
opts = f(opts, 'n_iter_admm', 10);
opts = f(opts, 'n_iter_cg', 15);
opts = f(opts, 'cost', {});
opts = f(opts, 'cg_tol', 1e-12);
opts = f(opts, 'cg_display_ind', []);


end