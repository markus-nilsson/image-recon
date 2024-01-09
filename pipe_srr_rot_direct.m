function data_srr = pipe_srr_rot_direct(data_lr, data_hr, lambda)
% function data_srr = pipe_srr_rot_direct(data_lr, data_hr, lambda)
%
% for data rotated around the AP axis (axis 2)

if (nargin < 3), lambda = 0.05; end

% compute data in hr space first

O1_list = op_obj_image_h2l.make_many(data_lr.h, data_hr.h, 3);
O1 = op_obj_image_stack(O1_list); 

data_tmp = O1' * data_lr;

% sharpen, start by defining this for one slice
tmp_lr = data_lr.copy(); tmp_lr.trim([], 1, [], []);
tmp_hr = data_hr.copy(); tmp_hr.trim([], 1, [], []);

O_list = op_obj_image_h2l.make_many(tmp_lr.h, tmp_hr.h, 3);
O = op_obj_image_stack(O_list); 

S = [];
for c = 1:numel(O_list)
    S = cat(1, S, O_list{c}.S);
end


NR = numel(O_list);
alpha = tmp_lr.data_obj{1}.pixdim(3) / tmp_lr.data_obj{1}.pixdim(1);

I = speye(size(S,2), size(S,2));
P = inv( (1-lambda) * (S' * S) + lambda * NR * alpha * I);


Q = data_tmp.imreshape();
for k = 1:data_tmp.dim(4)
    for j = 1:data_tmp.dim(2)
        T = squeeze(Q(:, j, :, k));
        T = P * T(:);
        Q(:, j, :, k) = reshape(T, data_tmp.dim(1), 1, data_tmp.dim(3));
    end
end

data_srr = data_hr.new(Q);


