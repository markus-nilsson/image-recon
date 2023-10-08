classdef op_obj_image_h2l < op_obj_image

    methods

        function O = op_obj_image_h2l(h_lr, h_hr, k)

            if (nargin < 3), k = 1; end

            % compute sampling matrix
            S = op_obj_image_h2l.S_h2l(h_hr, h_lr);

            O = O@op_obj_image(S, h_hr, h_lr, k);

            % somewhat temporary fix: address issue of sampling artifacts
            j = prod(h_lr.pixdim(2:5)) / prod(h_hr.pixdim(2:5));

            tmp = O.M * ones(size(O.M,2), 1);
            tmp = (1./tmp) * j;
            tmp(isinf(tmp)) = 0;
            tmp = spdiags(tmp, 0, numel(tmp), numel(tmp));
            
            O.aM = tmp * O.M;
            
            tmp = O.M' * ones(size(O.M,1), 1);
            tmp = (1./tmp) * 1/j;
            tmp(isinf(tmp)) = 0;
            tmp = spdiags(tmp, 0, numel(tmp), numel(tmp));
            
            O.aMT = tmp * O.M';

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
        
        function O_list = make_many(h_lrs, h_hr) % input: nifti headers
            O_list = cell(1, numel(h_lrs));
            for c = 1:numel(h_lrs)
                O_list{c} = op_obj_image_h2l(h_lrs{c}, h_hr);
            end
        end

    end
end
