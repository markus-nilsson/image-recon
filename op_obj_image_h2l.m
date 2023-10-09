classdef op_obj_image_h2l < op_obj_image

    methods

        function O = op_obj_image_h2l(h_lr, h_hr, n_j)

            % h_lr - header to low resolution image
            % h_hr - header to high resolution image
            % l    - index to operate on (w_jl, j - voxel index)

            % compute sampling matrix
            S = op_obj_image_h2l.S_h2l(h_hr, h_lr);

            O = O@op_obj_image(S, h_hr, h_lr, n_j);

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
                op_obj_image.h2ind(h_hr), h_hr), h_lr);
            
            S = op_obj_image.linear_interpolation(...
                ind_lr, h_hr.dim(2:4), h_lr.dim(2:4));
        end        
        
        function O_list = make_many(h_lrs, h_hr, n_k) % input: nifti headers
            O_list = cell(1, numel(h_lrs));
            for c = 1:numel(h_lrs)
                O_list{c} = op_obj_image_h2l(h_lrs{c}, h_hr, n_k);
            end
        end

    end
end
