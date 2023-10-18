classdef op_obj_image_distort < op_obj_image

    methods

        function O = op_obj_image_distort(t, h, n_k)

            % not sure where this code should go
            % ideally, null along a phase encode axis instead
            t = reshape(double(t), prod(size(t,1,2,3)), size(t,4));
%             t(:, 1) = 0;
%             t(:, 3) = 0;

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


        function [DF, I] = nonlinear_registration(data_a, data_b)

            % 1. Logistics
            op = msf_tmp_path(1);

            a_fn = fullfile(op, 'tmp_a.nii.gz');
            b_fn = fullfile(op, 'tmp_b.nii.gz');

            mdm_nii_write(data_a.imreshape(), a_fn, data_a.h);
            mdm_nii_write(data_b.imreshape(), b_fn, data_b.h);

            % 2. Registration
            p = elastix_p_affine(800);
            p.Transform = 'BSplineTransform';
            p.FinalGridSpacingInVoxels  = [ 1 1 1 ] * 20;
            p_fn = elastix_p_write(p, 'p.txt');

            opt.mio.coreg.clear_header = 0;
            opt.do_overwrite = 1;

            [o_fn,~,T] = mdm_coreg(a_fn, b_fn, p_fn, op, opt);

            I = mdm_nii_read(o_fn);

            % 3. Logistics -> Save displacemnt fields
            DF = elastix_get_deformation_field(T.t);

            if (size(DF, 5) ~= 3), error('using an old framework? bad nii_read'); end

            DF = permute(DF, [1 2 3 5 4]);

        end
    end
end