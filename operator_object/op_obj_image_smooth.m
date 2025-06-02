classdef op_obj_image_smooth < op_obj_image

    properties

        fwhm;

    end

    methods

        % reweigh the elements by .* operation

        function O = op_obj_image_smooth(fwhm, h_rhs, h_lhs, n_k)

            O = O@op_obj_image([], h_rhs, h_lhs, n_k);

            O.fwhm = fwhm;

            O.n_i = prod(h_lhs.dim(2:4));
            O.n_j = prod(h_rhs.dim(2:4));

        end

        function w = i_apply(O, d, ind)
            if (ind ~= 1), error('not tested'); end
            tmp = d.imreshape();
            tmp = mio_smooth_4d(tmp, O.fwhm);
            w = reshape(tmp, size(d.w));
        end

        function w = i_apply_adjoint(O, d, ind)
            if (ind ~= 1), error('not tested'); end
            tmp = d.imreshape();
            tmp = mio_smooth_4d(tmp, O.fwhm);
            w = reshape(tmp, size(d.w));
        end

    end

end