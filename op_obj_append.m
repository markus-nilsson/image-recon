classdef op_obj_append < op_obj

    properties

        A, B;

    end

    properties (Access = private)
        c_type;
    end
    
    methods
        
        function O = op_obj_append(A,B)
            
            % assume A and B are of types kernel or image sampler

            type_a = 'op_obj_image';
            type_b = 'op_obj_kernel';

            if (my_isa(B, 'op_obj_append'))
                B_tmp = B.B;
            else
                B_tmp = B;
            end

            if (my_isa(A, 'op_obj_append'))
                A_tmp = A.A;
            else
                A_tmp = A;
            end


            if (0) 
                1;
            elseif (my_isa(A_tmp, type_a) && my_isa(B_tmp, type_a))
                O.c_type = 1;
            elseif (my_isa(A_tmp, type_b) && my_isa(B_tmp, type_b))
                O.c_type = 2;
            elseif (my_isa(A_tmp, type_a) && my_isa(B_tmp, type_b)) 
                O.c_type = 3;
            elseif (my_isa(A_tmp, type_b) && my_isa(B_tmp, type_a)) 
                O.c_type = 4;
            else 
                error('check this');
            end

            switch (O.c_type)
                case 1 % two image samplers
%                     assert(A.n_k == B_tmp.n_k, 'n_k differs');
%                     assert(A.n_l == B_tmp.n_l, 'n_l differs');
                    O.n_i = B_tmp.n_i;
                    O.n_k = A.n_k;
                    O.n_j = A.n_j;
                    O.n_l = A.n_l; 
                case 2 % two dictionaries – no use case?
                    error('check this');
                case 3 % image sampling first, then kernel
                    
                    O.n_i = A.n_i;
                    O.n_j = A.n_j;

                    O.n_k = B.n_k;
                    O.n_l = B.n_l; 
                
                    
                case 4 % A: kernel first, B: image sampling
                    O.n_i = B_tmp.n_i;
                    O.n_j = B_tmp.n_j;
                    O.n_k = A.n_k;
                    O.n_l = A.n_l;
            end

            O.A = A;
            O.B = B;

        end

        function x = init_x(O, a, b)
            if (nargin < 2), a = O.n_j; end
            if (nargin < 3), b = O.n_k; end

            x = O.A.init_x(a, b);
            
            % if (~ ((my_isa(O.B, 'op_obj_image') || my_isa(O.B, 'op_obj_append'))))
            %     x = O.A.init_x(a, b);
            % else
            %     x = O.B.init_x(a, b);
            % end

        end
        
        function y = apply(O, x, ind)
            if (nargin < 3), ind = []; end
            tmp = O.A.apply(x, ind);
            y = O.B.apply(tmp);
        end

        function y = apply_adjoint(O, x, ind)
            if (nargin < 3), ind = []; end
            y = O.A.apply_adjoint(O.B.apply_adjoint(x, ind));
        end

        
    end
end
