classdef op_obj_image_stack < op_obj_image

    properties
        O_list = {};
    end

    methods

        function O = op_obj_image_stack(O_list)

            O = O@op_obj_image();

            O.O_list = O_list;

            O.n_vox = O.O_list{1}.n_vox;
            O.n_mp  = max(cellfun(@(x) max(x.k), O.O_list));
            
        end
       
        function data_hr = init_x(O)
            data_hr = O.O_list{1}.init_x();
        end

        % example use: make multiple samles of high resolution object
        % with different sampling orientations
        function data_lhs = apply(O, data_rhs)
            data_lhs = data_obj_cell(numel(O.O_list));
            for c = 1:numel(O.O_list)
                data_lhs.data_obj{c} = O.O_list{c}.apply(data_rhs);
            end

        end

        function data_rhs = apply_adjoint(O, data_lhs, ind)
            if (nargin < 3), ind = 1:numel(O.O_list); end

            assert(my_isa(data_lhs, 'data_obj_cell'), 'wrong type'); 

            data_rhs = O.init_x(); 
            for c = ind                
                S = O.O_list{c};
                data_rhs.add( S.apply_adjoint(data_lhs.data_obj{c}), S.k);
            end         
        end
    end
end
