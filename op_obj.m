classdef op_obj

    properties
        transpose = 0;
        n_i = 0; % low-res size
        n_j = 0; % high-res size (fixed)
        n_k = 0; % model coefficients (fixed)
        n_l = 0; % experimental or imaging settings

        h_lhs;
        h_rhs;

    end
    
    methods
        
        function O = op_obj()
            O.h_lhs = [];
            O.h_rhs = [];
        end
        
        function O = ctranspose(O)
            O.transpose = ~O.transpose;
        end
        
        function y = mtimes(O,x)

            assert(my_isa(O, 'op_obj'), 'lhs must be an op_obj');

            if (~O.transpose)
                y = O.apply(x);
            else
                y = O.apply_adjoint(x);
            end
            
        end

        % this should init x with zeros; smarter initalizations may be
        % possible, but implementing this within this function would break
        % things (a better function name would have been good, but hey...)
        function data_rhs = init_x(O, a, b) %
            if (nargin < 2), a = O.n_j; end
            if (nargin < 3), b = O.n_k; end

            if (~isempty(O.h_rhs))
                h = O.h_rhs;
                h.dim(5) = b;
                data_rhs = do_w_image_vector(zeros(a,b), h);
            else
                data_rhs = do_w(zeros(a,b));
            end
        end
        

    end

    methods (Abstract)
        y = apply(O,x);
        y = apply_adjoint(O,x);
    end


end
