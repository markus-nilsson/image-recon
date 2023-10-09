classdef op_obj

    properties
        transpose = 0;
        n_i = 0; % low-res size
        n_j = 0; % high-res size (fixed)
        n_k = 0; % model coefficients (fixed)
        n_l = 0; % experimental or imaging settings
    end
    
    methods
        
        function O = op_obj()
            if (nargin == 0), return; end
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

    end

    methods (Abstract)
        y = apply(O,x);
        y = apply_adjoint(O,x);
        x = init_x(O);
    end


end
