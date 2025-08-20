function data_ref = mec_ref_b0(data, opts)

if (nargin < 2), opts.present = 1; end

data_ref = data.copy();

data_ref.w = repmat(...
    mean(data_ref.w(:,data.xps.b < 0.05e9),2), 1, data.xps.n);

sc = diag(data.w' * data_ref.w) ./ diag(data.w' * data.w);

data_ref.w = data_ref.w * diag(1./sc);
