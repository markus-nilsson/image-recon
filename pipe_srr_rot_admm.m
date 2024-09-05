function [data_srr, O] = pipe_srr_rot_admm(data_lr, data_hr, lambda)
% function data_srr = pipe_srr_rot_direct(data_lr, data_hr, lambda)
%
% for data rotated around the AP axis (axis 2)

% compute data in hr space first

if (nargin < 3), lambda = 0.05; end

if (~all(data_lr.dim(4) == median(data_lr.dim(4))))
    error('check data, dim(4) varies');
end


n_k = data_lr.dim(4);
if (~isequal(n_k, n_k-n_k + n_k(1)))
    error('check this');
end
n_k = n_k(1);

O_list = op_obj_image_h2l.make_many(data_lr.h, data_hr.h, n_k);
O = op_obj_image_stack(O_list);

ind_flt = ones(1, min(data_lr.dim(4)));

opts.cost = {...
    cost_imfilter_3d(1, lambda, ind_flt), ...
    cost_positivity(1, ind_flt), ...
    };
opts.n_iter_admm = 5;
opts.n_iter_cg = 20;

[data_srr, r, cgr] = srr_est_admm(O, data_lr, opts);

1;
