classdef do_c < do

    properties
        data_obj = {};   
        h;
    end

    methods

        % data object as cell - the data object w_il is here represented
        % as (w_l){l}, meaning the second dimension is the number of cells
        % 
        % note that i = 1..n_vox, where n_vox(l)        

        function O = do_c(varargin)

            if (nargin == 0), return; end

            a = varargin{1};
            if (isnumeric(a))
                O.data_obj = cell([1 a]);
            elseif (iscell(a)) && (my_isa(a{1}, 'do_w'))
                O.data_obj = a;
            end

        end

        function O_new = new(O, data_obj)
            O_new = do_c(data_obj);
        end

        function O_new = copy(O)
            tmp = cell(size(O.data_obj));
            for c = 1:numel(O.data_obj)
                tmp{c} = O.data_obj{c}.copy();
            end
            O_new = new(O, tmp);
        end

        function O_new = subsample(O, ind)
            O_new = do_c(O.data_obj(ind));
        end

        function O = select(O, ind)
            if (numel(ind)~=1), error('assume selecting single object'); end
            O = O.data_obj{ind};
        end

        function I = imreshape(O)

            A = zeros(4, numel(O.data_obj));
            for c = 1:numel(O.data_obj)
                A(:,c) = O.data_obj{c}.dim(1:4);
            end

            if (~all(A == A(:,1)))
                error('different dims of data objs');
            end

            if (~all(A(4,:) == 1))
                error('cannot deal with 4d objs yet');
            end

            I = zeros(A(:,1)');

            for c = 1:numel(O.data_obj)
                I(:,:,:,c) = O.data_obj{c}.imreshape();
            end


        end
        

        function q = minus(a,b)

            ac = isa(a,'do_c');
            bc = isa(b,'do_c');

            na = numel(a);
            nb = numel(b);

            if (~ac) || (~bc), error('expected both to be data_obj_cell'); end
            if (na ~= nb), error('different size of the two objs'); end
            
            q = do_c(max(na, nb));
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

        function O = trim(O, ir, jr, kr, vr)
            for c = 1:numel(O)
                O.data_obj{c}.trim(ir, jr, kr, vr);
            end
        end

        function dim = dim(O, c_dim)
            if (nargin < 2), c_dim = 1:4; end
            dim = zeros(numel(c_dim), numel(O.data_obj));
            for c = 1:size(dim,2)
                dim(:, c) = O.data_obj{c}.dim(c_dim);
            end
        end
        
        function h = get.h(O)
            h = cell(1, numel(O));
            for c = 1:numel(h)
                h{c} = O.data_obj{c}.h;
            end
        end

        

        function s = sum(O)
            error('not implemented');
        end



    end
end

