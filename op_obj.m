classdef op_obj

    properties
        M;
        n_vox;
        n_mp;
        k;

        transpose = 0;        
    end
    
    methods
        
        function O = op_obj(M)
            if (nargin == 0), return; end
            O.M = M;
            O.n_vox = size(M,2);
            O.n_mp = 1;
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

        function x = init_x(O)
            x = data_obj(zeros(O.n_mp, O.n_vox));
        end

        function y = apply(O,x)
            y = O.M * op_obj.f(x);
        end

        function y = apply_adjoint(O,x)
            y = O.M' * op_obj.f(x);
        end

    end

end
