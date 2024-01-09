classdef op_obj_image_logexp < op_obj_image

    properties

        S;

    end

    methods

        % reweigh the elements by .* operation

        function O = op_obj_image_logexp(S, h_rhs, h_lhs, n_k)

            O = O@op_obj_image([], h_rhs, h_lhs, n_k);

            O.S = S.^2;

        end

        function x = init_x(O, a, b)
            if (nargin < 2), a = O.n_j; end
            if (nargin < 3), b = O.n_k; end
            h = O.h_rhs;
            h.dim(5) = b;
            x = do_w_image_log(zeros(a,b), h);
        end

        function w = i_apply(O, d, ind)
            w = exp(d.w(:, ind));
        end

        function w = i_apply_adjoint(O, d, ind)
            w = real(log(d.w(:, ind)));
            w(isinf(w(:))) = 0;
        end


    end

end