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

        function w = i_apply(O, d, ind)
            if (size(O.W,2) == 1) && (max(ind) > 1) % temp fix
                w = O.W .* d.w(:, ind);
            else
                w = O.W(:, ind) .* d.w(:, ind);
            end

        end

        function w = i_apply_adjoint(O, d, ind)
            if (size(O.WT,2) == 1) && (max(ind) > 1) % temp fix
                w = O.WT .* d.w(:, ind);
            else
                w = O.WT(:, ind) .* d.w(:, ind);
            end

        end

    end

end