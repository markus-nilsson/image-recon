classdef op_obj_image_distort < op_obj_image

    methods

        function O = op_obj_image_distort(t, h, n_k)

            % not sure where this code should go
            % ideally, null along a phase encode axis instead
            t = reshape(double(t), prod(size(t,1,2,3)), size(t,4));
            t(:, 1) = 0;
            t(:, 3) = 0;

            % compute sampling matrix
            S = op_obj_image_distort.S_distort(t, h);

            O = O@op_obj_image(S, h, h, n_k);

        end

    end

    methods (Static)

        function S = S_distort(t, h)

            t = cat(2, t, zeros(size(t,1), 1))';

            pc = op_obj_image.ind2pc(...
                op_obj_image.h2ind(h),h);
            pc = pc + t;

            ind = op_obj_image.pc2ind(pc, h);

            dim = h.dim(2:4)';

            S = op_obj_image.linear_interpolation(ind, dim, dim);
        end        

    end
end
