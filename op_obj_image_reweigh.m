classdef op_obj_image_reweigh < op_obj_image

    properties

        W;
        WT;

    end

    methods

        % reweigh the elements by .* operation

        function O = op_obj_image_reweigh(W, WT, h_rhs, h_lhs, n_k)

            O = O@op_obj_image([], h_rhs, h_lhs, n_k);

            O.n_i = prod(h_lhs.dim(2:4));
            O.n_j = prod(h_rhs.dim(2:4));

            O.W = W;
            O.WT = WT;


        end

        function data_lhs = apply(O, data_rhs, ind)
            if (nargin < 3) || (isempty(ind)), ind = 1:size(data_rhs.w,2); end
            assert(my_isa(data_rhs, O.d_type), 'rhs not a %s', O.d_type);
            data_lhs = do_w_image_vector(O.W(:,ind) .* data_rhs.w(:,ind), O.h_lhs);
        end

            


        function data_rhs = apply_adjoint(O, data_lhs, ind)
            if (nargin < 3) || (isempty(ind)), ind = 1:size(data_lhs.w,2); end
            assert(my_isa(data_lhs, O.d_type), 'rhs not a %s', O.d_type);
            data_rhs = do_w_image_vector(O.WT(:,ind) .* data_lhs.w(:,ind), O.h_rhs);
        end

    end

end