classdef cost_imprior < cost_admm

    properties

        c_type;
        z;

        im;

    end

    methods

        function C = cost_imprior(mu, im, ind)
            if (nargin < 2), ind = []; end
            C = C@cost_admm(mu, ind);
            C.im = im;
        end

        function [f,b] = do_iter(C, x)

            C.z = C.image_filter(x, C.ind);

            [f,b] = C.admm_update(x, C.z);

        end

        function x_flt = image_filter(C, x, ind)

            if (all(all(x.w == 0)))
                x_flt = x.new(x.w);
                return;
            end

            P = my_clean_prediction(C.im.imreshape(), x.imreshape(), 0);      

            x_flt = x.new( P );

        end


    end


    methods (Static)

        function p = morph_prior(m,y)

            d = -4:4;

            for i = 1:9
                for j = 1:9
                    for k = 1:9
                        di(i,j,k) = (i-5);
                        dj(i,j,k) = (j-5);
                        dk(i,j,k) = (k-5);
                    end
                end
            end

            p = m.zeros();

            for i = 1:m.dim(1)
                for j = 1:m.dim(2)
                    for k = 1:m.dim(3)

                        ir = i + di;
                        jr = j + dj;
                        kr = k + dj;

                        mind = (ir >= 1) & (jr >= 1) & (kr >= 1) & ...
                            (ir <= m.dim(1)) & ...
                            (jr <= m.dim(2)) & ...
                            (kr <= m.dim(3));

                        ir = ir(mind);
                        jr = jr(mind);
                        kr = kr(mind);

                        ind0 = sub2ind(m.dim(1:3), i, j, k);
                        ind  = sub2ind(m.dim(1:3), ir, jr, kr);

                        a = m.w(ind);
                        b = y.w(ind, 1);

                        w = (a - m.w(ind0)).^2;
                        w = exp(-w.^2 / 25^2);
                        w = w / sum(w(:));

                        p.w(ind0) = median( b(w(:) > median(w(:))));
                    end
                end
            end
        end
    end
end