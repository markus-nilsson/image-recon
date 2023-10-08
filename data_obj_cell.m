classdef data_obj_cell < handle

    properties
        data_obj = {};
    end

    methods

        function O = data_obj_cell(varargin)

            if (nargin == 0), return; end

            a = varargin{1};
            if (isnumeric(a))
                O.data_obj = cell([1 a]);
            elseif (iscell(a)) && (my_isa(a{1}, 'data_obj'))
                O.data_obj = a;
            end

        end

        function O_new = subsample(O, ind)
            O_new = data_obj_cell(O.data_obj(ind));
        end
        

        function q = minus(a,b)

            ac = isa(a,'data_obj_cell');
            bc = isa(b,'data_obj_cell');

            na = numel(a);
            nb = numel(b);

            if (~ac) || (~bc), error('expected both to be data_obj_cell'); end
            if (na ~= nb), error('different size of the two objs'); end
            
            q = data_obj_cell(max(na, nb));
            for c = 1:numel(q)
                q.data_obj{c} = a.data_obj{c} - b.data_obj{c};
            end

        end

        function n = numel(O)
            n = numel(O.data_obj);
        end

        function r = norm(O)
            r = 0;
            for c = 1:numel(O.data_obj)
                r = r + norm(O.data_obj{c}); 
            end
        end

        function O = add_noise_gaussian(O, noise_std)
            for c = 1:numel(O.data_obj)
                O.data_obj{c}.add_noise_gaussian(noise_std);
            end
        end

    end
end

