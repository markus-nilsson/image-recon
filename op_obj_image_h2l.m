classdef op_obj_image_h2l < op_obj_image

    methods

        function O = op_obj_image_h2l(h_lr, h_hr, n_k)

            % h_lr - header to low resolution image
            % h_hr - header to high resolution image
            % l    - index to operate on (w_jl, j - voxel index)

            % compute sampling matrix
            S = op_obj_image_h2l.S_h2l(h_hr, h_lr);

            O = O@op_obj_image(S, h_hr, h_lr, n_k);

            % somewhat temporary fix: address issue of sampling artifacts
            j = prod(h_lr.pixdim(2:5)) / prod(h_hr.pixdim(2:5));

            tmp = S * ones(size(S,2), 1);
            tmp = (1./tmp) * j;
            tmp(isinf(tmp)) = 0;
            tmp = spdiags(tmp, 0, numel(tmp), numel(tmp));

            O.aM = tmp * S;

            tmp = S' * ones(size(S,1), 1);
            tmp = (1./tmp) * 1/j;
            tmp(isinf(tmp)) = 0;
            tmp = spdiags(tmp, 0, numel(tmp), numel(tmp));

            O.aMT = tmp * S';

        end

    end

    methods (Static)

        function S = S_h2l(h_hr, h_lr)
            % Compute the high-to-low resolution operator

            ind_lr = op_obj_image.pc2ind(...
                op_obj_image.ind2pc(...
                op_obj_image.h2ic(h_hr), h_hr), h_lr);

            S = op_obj_image.linear_interpolation(...
                ind_lr, h_hr.dim(2:4), h_lr.dim(2:4));
        end

        function O_list = make_many(h_lrs, h_hr, n_k) % input: nifti headers
            O_list = cell(1, numel(h_lrs));
            for c = 1:numel(h_lrs)
                O_list{c} = op_obj_image_h2l(h_lrs{c}, h_hr, n_k);
            end
        end


        function S = tmp_rotate2d(sz_hr, theta, ar)
            % temporary function, should be rewritten (now used in examples
            % only)

            if (numel(ar) ~= 1), error('one ar only'); end

            % downsampling in y direction and rotation
            R_h2l = diag([1 1/ar 1]);
            R_rot = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

            R_tr_hr  = [1 0 -sz_hr(1)/2; 0 1 -sz_hr(2)/2; 0 0 1];
            R_tr_lr  = [1 0 +sz_hr(1)/2; 0 1 +sz_hr(2)/2/ar; 0 0 1];

            R = R_h2l * inv(R_tr_hr) * R_rot * R_tr_hr;

            sz_lr = ceil(sz_hr ./ [1 ar]);


            S = spalloc(prod(sz_lr), prod(sz_hr), prod(sz_hr)*4);

            for i_hr = 1:sz_hr(1)
                for j_hr = 1:sz_hr(2)

                    p = R * [i_hr j_hr 1]';

                    p0 = floor(p);
                    w1 = p(1) - p0(1);
                    w2 = p(2) - p0(2);

                    w = [ (1-w1)*(1-w2)  w1*(1-w2)  (1-w1)*w2  w1*w2];

                    p = p0(1:2)' + [0 0; 1 0; 0 1; 1 1];

                    ind = ...
                        (p(:,1) > 0) & (p(:,1) <= sz_lr(1)) & ...
                        (p(:,2) > 0) & (p(:,2) <= sz_lr(2));

                    i_lr = p(ind,1);
                    j_lr = p(ind,2);

                    S( (j_lr-1)*sz_lr(1) + i_lr, (j_hr-1)*sz_hr(1) + i_hr) = w(ind);

                end
            end


        end
    end
end