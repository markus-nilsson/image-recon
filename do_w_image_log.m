classdef do_w_image_log < do_w_image_vector

    properties 
        S;
    end

    methods

        function O = do_w_image_log(w, h, do_log)

            if (nargin < 3), do_log = 1; end

            if (do_log)

                S = w;
                S(S < eps) = eps;

                w = log(S);

                S(S < 10  * eps) = 0;

                S = S.^2;
                S = S / mean(S(:));

            else
                S = [];
            end
            
            O@do_w_image_vector(w, h);
            O.S = S;

        end

        function O_new = new(O, w)
            O_new = do_w_image_log(w, O.h, 0);
            O_new.S = O.S;
        end

        function r = norm(O)
            r = O.w.^2;
%             r = r .* O.S;
            r(isnan(r(:))) = 0;
            r = sum(r,2);
        end


    end

end

