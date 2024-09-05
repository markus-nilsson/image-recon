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

            O = O@do();
            O.c_type = 2;

            if (nargin == 0), return; end

            a = varargin{1};
            if (isnumeric(a))
                O.data_obj = cell([1 a]);
            elseif (iscell(a)) && (my_isa(a{1}, 'do_w'))
                O.data_obj = a;
            end

        end

        function O_new = new(O, data_obj)
            
            if (iscell(data_obj))
                O_new = do_c(data_obj);
                O_new.h = O.h;
            elseif (isnumeric(data_obj))
                O_new = do_c(numel(O));
                O_new.h = O.h;
                ind0 = 0;
                for c = 1:numel(O_new)
                    n = size(O.data_obj{c}.w,2);
                    ind = (1:n) + ind0;
                    ind0 = ind0 + n;
                    switch (ndims(data_obj))
                        case 2
                            TMP = data_obj(:,ind);
                        case 4
                            TMP = data_obj(:,:,:,ind);
                        otherwise
                            error('check your data');
                    end
                    O_new.data_obj{c} = O.data_obj{c}.new(TMP);
                end
            end
        end

        % slower implementation than necessary (creating big matrix twice)
        function O_new = zeros(O)
            O_new = O.new(zeros(size(O.imreshape())));
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

        % warning: not implemented yet
        function O_new = conj(O)
            O_new = O; 
        end

        function O = select(O, ind)
            if (numel(ind)~=1), error('assume selecting single object'); end
            O = O.data_obj{ind};
        end

        function w = flatten(O)
            TMP = O.imreshape();
            w = reshape(TMP, prod(size(TMP,1,2,3)), size(TMP,4));
        end

        function I = imreshape(O)

            A = zeros(4, numel(O.data_obj));
            for c = 1:numel(O.data_obj)
                A(:,c) = O.data_obj{c}.dim(1:4);
            end

            if (~all(A(1,:) == A(1,1))) || ...
                    (~all(A(2,:) == A(2,1))) || ...
                    (~all(A(3,:) == A(3,1)))
                error('different dims of data objs');
            end

            I = zeros(A(1,1), A(2,1), A(3,1), sum(A(4,:)));

            ind0 = 0;
            for c = 1:numel(O.data_obj)
                TMP = O.data_obj{c}.imreshape();
                ind = (1:size(TMP,4)) + ind0;
                ind0 = ind0 + size(TMP,4);
                I(:,:,:,ind) = TMP;
            end

        end
        

        function q = mtimes(a,b)
            q = apply_operation(a, b, 1);
        end

        function q = plus(a,b)
            q = apply_operation(a, b, 2);
        end
        
        function q = minus(a,b)
            q = apply_operation(a, b, 3);
        end

        function q = times(a,b)
            q = apply_operation(a, b, 4);
        end

        function q = apply_operation(a,b,c_type)

            if (isnumeric(a))
                fa = @(c) a;
            elseif (isa(a, 'do_c'))
                fa = @(c) a.data_obj{c};
            else
                error('check this');
            end

            if (isnumeric(b))
                fb = @(c) a;
            elseif (isa(b, 'do_c'))
                fb = @(c) b.data_obj{c};
            else
                error('check this');
            end
            
            q = do_c(max(numel(a), numel(b)));
            for c = 1:numel(q)
                switch  (c_type)
                    case 1
                        q.data_obj{c} = fa(c) * fb(c);
                    case 2
                        q.data_obj{c} = fa(c) + fb(c);
                    case 3
                        q.data_obj{c} = fa(c) - fb(c);
                    case 4
                        q.data_obj{c} = fa(c) .* fb(c);
                end
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
            s = 0;
            for c = 1:numel(O)
                s = s + sum(O.data_obj{c});
            end
        end

    end
end

