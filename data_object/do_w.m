classdef do_w < do

    properties
        n_vox = 0;
        w = [];
        xps = [];

    end

    methods

        % do_w is a data object with weights, given as a w_il matrix
        % where i of the numbe of voxels and l the number of model
        % coefficients or contrasts

        function O = do_w(w)

            O = O@do();
            O.c_type = 1;

            if (canUseGPU)
                O.w = gpuArray(single(w)); 
            else
                O.w = double(w); % necessary for sparse operations
            end
            
            O.n_vox = size(O.w, 1);
        end

        function O_new = new(O, w)
            O_new = do_w(w);
            O_new.transpose = O.transpose;
        end

        function O_new = copy(O)
            O_new = O.new(O.w);
        end

        function O_new = subsample(O, ind)
            O_new = O.new(O.w(:, ind));
        end

        function O_new = mean(O, ind)
            O_new = O.new(mean(O.w(:, ind),2));
        end

        function O_new = min(O, ind)
            O_new = O.new(min(O.w(:, ind),[],2));
        end

        function q = minus(a,b)
            q = do_w.apply_operator(a,b,@minus);
        end

        function q = plus(a,b)
            q = do_w.apply_operator(a,b,@plus);
        end

        function q = times(a,b)
            q = do_w.apply_operator(a,b,@times);
        end

        function q = mtimes(a,b)
            q = do_w.apply_operator(a,b,@mtimes);
        end

        function r = norm(O)
            r = sum(O.w(~isnan(O.w(:))).^2);
        end

        function s = sum(O)
            s = sum(O.w(:));
        end

        function n = numel(O)
            n = numel(O.w(:));
        end

        function N = conj(O)
            N = O.new(real(O.w) - 1i * imag(O.w));
        end

        function O = add_gaussian_noise(O, noise_std)
            f = @(x) x + randn(size(x)) * noise_std;
            O.w = f(O.w);
        end

    end

    methods (Static)

        function q = apply_operator(a,b,op)

            ca = my_isa(a, 'do_w');
            cb = my_isa(b, 'do_w');

            if (ca && cb)
                q = a.new(op(a.w, b.w));
            elseif (ca) && (isnumeric(b))
                q = a.new(op(a.w, b));
            elseif (cb) && (isnumeric(a))
                q = b.new(op(a, b.w));
            else
                error('not implemented');
            end

        end

    end

end

