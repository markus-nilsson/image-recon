classdef data_obj < handle

    properties
        n_vox = 0;
        w = [];

        transpose = 0;
    end

    methods

        function O = data_obj(w)
            O.w = double(w); % necessary for sparse operations
            O.n_vox = size(O.w, 1);
        end

        function O_new = new(O, w)
            O_new = data_obj(w);
            O_new.transpose = O.transpose;
        end

        function O = ctranspose(O)
            O.transpose = ~O.transpose;
        end

        function q = minus(a,b)
            q = data_obj.apply_operator(a,b,@minus);
        end

        function q = plus(a,b)
            q = data_obj.apply_operator(a,b,@plus);
        end

        function q = times(a,b)
            q = data_obj.apply_operator(a,b,@times);
        end

        function q = mtimes(a,b)
            q = data_obj.apply_operator(a,b,@mtimes);
        end

        function r = norm(O)
            r = sum(O.w(:).^2);
        end

        function s = sum(O)
            s = sum(O.w(:));
        end

        function n = numel(O)
            n = numel(O.w(:));
        end

        function O = add_gaussian_noise(O, noise_std)
            f = @(x) x + randn(size(x)) * noise_std;
            O.w = f(O.w);
        end        
        

    end

    methods (Static)

        function q = apply_operator(a,b,op)

            ca = my_isa(a, 'data_obj');
            cb = my_isa(b, 'data_obj');

            if (ca && cb)
                q = a.new(op(a.w, b.w)); 
            elseif (ca)
                q = a.new(op(a.w, b)); 
            elseif (cb)
                q = b.new(op(a, b.w)); 
            else
                error('unexpected');
            end

        end

    end

end

