classdef op_obj_image_stack < op_obj_image

    properties
        O_list = {};
        l_list = [];
    end

    methods

        % acts on data objects as cells (do_c)

        function O = op_obj_image_stack(O_list, l_list)

            if (nargin < 2)
                for c = 1:numel(O_list)
                    l_list{c} = 1:O_list{c}.n_l;
                end
            end

            assert(all(O_list{1}.n_j == cellfun(@(x) x.n_j, O_list)), ...
                'All objects must act on the same hr object (n_j)');            

            assert(all(O_list{1}.n_k == cellfun(@(x) x.n_k, O_list)), ...
                'All objects must act on the same hr object (n_k)');            
            
            O = O@op_obj_image();

            O.O_list = O_list;
            O.l_list = l_list;

            O.n_i = cellfun(@(x) x.n_i, O_list, 'UniformOutput', false);
            O.n_j = O.O_list{1}.n_j;
            O.n_k = O.O_list{1}.n_k;

            if (all(cellfun(@(x) isempty(x), O.l_list)))
                % temporary fix
                O.n_l = O.O_list{1}.n_l;
            else
                O.n_l = max(cellfun(@(x) max(x), O.l_list));
            end
            
        end
       
        function data_hr = init_x(O, a, b)
            if (nargin < 2), a = O.n_j; end
            if (nargin < 3), b = O.n_l; end 
            data_hr = O.O_list{1}.init_x(a, b);
        end

        % example use: make multiple samles of high resolution object
        % with different sampling orientations (e.g. many low-res objects)
        function data_lhs = apply(O, data_rhs, ind)

            if (nargin < 3), ind = []; end
            
            if (isempty(ind))
                ind = cell2mat(O.l_list')';
            end

            if (isempty(ind))
                f_ind = @(c) [];
            else
                f_ind = @(c) ind(:,c);
            end
            
            data_lhs = do_c(numel(O.O_list));
            for c = 1:numel(O.O_list)

                if (isa(data_rhs, 'do_c'))
                    % trying this, not the same behaviour, but one we want
                    data_lhs.data_obj{c} = O.O_list{c}.apply(data_rhs.data_obj{c}, f_ind(c));
                else                
                    data_lhs.data_obj{c} = O.O_list{c}.apply(data_rhs, f_ind(c));
                end


                % data_lhs.data_obj{c} = O.O_list{c}.apply(data_rhs, f_ind(c));
            end
        end

        function data_rhs = apply_adjoint(O, data_lhs, ind)
            
            if (nargin < 3) || (isempty(ind))
                ind = 1:numel(O.O_list); 
            end

            assert(my_isa(data_lhs, 'do_c'), 'wrong type'); 

            data_rhs = O.init_x(); 
            for c = ind          
                tmp = O.O_list{c}.apply_adjoint(data_lhs.data_obj{c});
                data_rhs.add( tmp, O.l_list{c});
            end         
        end
    end
end
