function mp = smr_admm(O, data, opts)

if (nargin < 3), opts.present = 1; end

opts = msf_ensure_field(opts, 'n_iter_admm', 40);
opts = msf_ensure_field(opts, 'n_iter_cg', 40);
opts = msf_ensure_field(opts, 'tol', 0);
opts = msf_ensure_field(opts, 'mu1', 1); 


O.mp = O.reset_mp();

O.c_type = 1;
C = cost_dict(O, opts.mu1);
C.data = data;
opts.cost_obj = {C};
[mp1,r1,cgr1] = smr_est_admm(O, data, opts);
O.mp(O.k02_ind) = mp1;

O.c_type = 2;
opts.cost_obj = {};
[mp2,r2,cgr2] = smr_est_admm(O, data, opts);
O.mp(O.odf_ind) = mp2;

O.c_type = 1;
C = cost_dict(O, opts.mu1);
C.data = data;
opts.cost_obj = {C};
[mp3,r3,cgr3] = smr_est_admm(O, data, opts);
O.mp(O.k02_ind) = mp3;

O.c_type = 2;
opts.cost_obj = {};
[mp4,r4,cgr4] = smr_est_admm(O, data, opts);
O.mp(O.odf_ind) = mp4;

O.c_type = 1;
C = cost_dict(O, opts.mu1);
C.data = data;
opts.cost_obj = {C};
[mp5,r5,cgr5] = smr_est_admm(O, data, opts);
O.mp(O.k02_ind) = mp5;

mp(O.odf_ind) = mp4;
mp(O.k02_ind) = mp5;
