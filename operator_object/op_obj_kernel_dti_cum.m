classdef op_obj_kernel_dti_cum < op_obj_kernel 

    properties
        S;
        h;
    end
    
    methods
        
        function O = op_obj_kernel_dti_cum(data)

            M = [ones(data.xps.n, 1) -data.xps.bt * 1e-9]';

            O = O@op_obj_kernel(M);

            O.h = data.h;
            O.n_j = prod(data.h.dim(2:4));
            O.n_i = O.n_j;
            O.n_k = 7; 

        end   

        function x = init_x(O, a, b)
            if (nargin < 2), a = O.n_j; end
            if (nargin < 3), b = O.n_k; end
            x = do_w_image_log(zeros(a, b), O.h);
        end
    
    end

end
