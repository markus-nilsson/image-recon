classdef op_obj_image_distort < op_obj_image

    methods

        function O = op_obj_image_distort(d, n_k)

            assert(d.dim(4) == 3, 'expected a displaceement field (d.dim(4)==3)');

            % !!!!
            % this is not a validated function
            % !!!!

            % t is the displacement field
            t = d.imreshape();
            t = reshape(t, prod(size(t,1,2,3)), size(t,4));

            % compute sampling matrix
            S = op_obj_image_distort.S_distort(t, d.h);

            O = O@op_obj_image(S, d.h, d.h, n_k);

        end

    end

    methods (Static)

        function S = S_distort(t, h)

            t = cat(2, t, zeros(size(t,1), 1));

            pc = op_obj_image.ind2pc(...
                op_obj_image.h2ic(h),h);

            pc = pc + t;

            ind = op_obj_image.pc2ind(pc, h);

            dim = h.dim(2:4)';

            S = op_obj_image.linear_interpolation(ind, dim, dim);
        end

        function o = interpolate(c, dim, ss)

            c = reshape(c', [dim(1) dim(2) dim(3) 4]);

            o = zeros([floor(dim(:)/ss); 4]'); 
            for c_l = 1:4
                o(:,:,:,c_l) = imresize3(c(:,:,:,c_l), 1/ss);
            end

            o = reshape(o, [prod(size(o, [1 2 3])) 4])';

        end


        function [DF, I, FD] = nonlinear_registration(data_a, data_b)

            % 1. Logistics
            op = msf_tmp_path(1);

            a_fn = fullfile(op, 'tmp_a.nii.gz');
            b_fn = fullfile(op, 'tmp_b.nii.gz');

            mdm_nii_write(data_a.imreshape(), a_fn, data_a.h);
            mdm_nii_write(data_b.imreshape(), b_fn, data_b.h);

            % 2. Registration
            p = elastix_p_affine(1000);
            p.Transform = 'BSplineTransform';
            p.FinalGridSpacingInVoxels  = [ 20 20 20 ] * 2.5;% * 20;
            p_fn = elastix_p_write(p, 'p.txt');

            opt.mio.coreg.clear_header = 0;
            opt.do_overwrite = 1;
            
            [o_fn,~,T] = mdm_coreg(a_fn, b_fn, p_fn, op, opt);

            if (1)
                I = mdm_nii_read(o_fn);
                A = data_a.imreshape();
                B = data_b.imreshape();

                msf_imagesc(cat(1, A, B, I));
                pause(0.05);
            end

            % Get the inverse transform
            p.Metric = 'DisplacementMagnitudePenalty';
            p.HowToCombineTransforms = 'Compose';
            p_fn = elastix_p_write(p, 'p.txt');
            
            [~,~,TI] = mdm_coreg(a_fn, a_fn, p_fn, op, opt, T.t);


            % 3. Logistics -> Save displacemnt fields
            DF = elastix_get_deformation_field(T.t);

            TI.t.InitialTransform = 'NoInitialTransform';
            TI.t = msf_rmfield(TI.t, 'InitialTransformParametersFileName');
            FD = elastix_get_deformation_field(TI.t);

            if (size(DF, 5) ~= 3), error('using an old framework? bad nii_read'); end

            DF = permute(DF, [1 2 3 5 4]);
            FD = permute(FD, [1 2 3 5 4]);

        end
    end
end