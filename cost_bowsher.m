classdef cost_bowsher < cost_admm

    properties

        z;
        S;

    end

    methods

        function C = cost_bowsher(mu, data, ind)
            if (nargin < 2), ind = []; end
            ind = logical(ind);
            C = C@cost_admm(mu, ind);


            
            if (size(data.dim(4)) ~= 1), error('must be a 3d volue'); end


            dind = [...
                1 0 0;
                0 1 0;
                0 0 1;
                -1 0 0;
                0 -1 0;
                0 0 -1;
                1 1 0;
                1 0 1;
                0 1 1;
                -1 -1 0;
                -1 0 -1;
                0 -1 -1;
                1 -1 0;
                0 1 -1;
                -1 0 1;
                -1 1 0;
                0 -1 1;
                1 0 -1]; % 18 neighbours

            C.S = ones(size(data.w,1), 4);
            for i = 1:data.dim(1)
                for j = 1:data.dim(2)
                    for k = 1:data.dim(3)

                        mind = ...
                            (dind(:,1) >= (i-1)) & ...
                            (dind(:,2) >= (j-1)) & ...
                            (dind(:,3) >= (k-1)) & ...
                            (dind(:,1) <= (data.dim(1)-i)) & ...
                            (dind(:,2) <= (data.dim(2)-j)) & ...
                            (dind(:,3) <= (data.dim(3)-k));

                        if (sum(mind) < 4), continue; end

                        c = [i j k] + dind(mind,:);

                        ind0 = sub2ind(data.dim(1:3), i, j, k);
                        ind  = sub2ind(data.dim(1:3), c(:,1), c(:,2), c(:,3));

                        [~,ind1] = sort(abs(data.w(ind0) - data.w(ind)), 1, 'ascend');

                        C.S(ind0, :) = ind(ind1(1:4));

                    end
                end
            end


        end

        function [f,b] = do_iter(C, x)

            C.z = C.image_filter(x, C.ind);

            [f,b] = C.admm_update(x, C.z);

        end

        function x_flt = image_filter(C, x, ind)

            x_flt = x.new(zeros(size(x.w)));
            
            for c = 1:find(ind)
                for d = 1:4
                    x_flt.w(:,c) = x_flt.w(:,c) + x.w(C.S(:,d), c) / 4;
                end
            end
        end


    end
end