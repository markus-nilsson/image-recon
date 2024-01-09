classdef op_obj_image_ft < op_obj_image

    properties (Access = protected)
        
    end

    methods

        function O = op_obj_image_ft(h_rhs, h_lhs, n_k)
            O = O@op_obj_image([], h_rhs, h_lhs, n_k);          
        end

        function w = i_apply(O, d, ind)
            w = d.imreshape();
            w = w(:,:,:,ind);
            w = fftn(w);
            w = reshape(w, prod(size(w,1,2,3)), size(w,4));
        end


        function w = i_apply_adjoint(O, d, ind)
            w = d.imreshape();
            w = w(:,:,:,ind);
            w = ifftn(w);
            w = reshape(w, prod(size(w,1,2,3)), size(w,4));
        end

    end


    methods (Static)

    end

end
