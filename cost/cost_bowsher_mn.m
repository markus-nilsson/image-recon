classdef cost_bowsher_mn < cost_admm

    % implement something bowsher-like
    properties

        z;
        S;

    end

    methods

        function C = cost_bowsher_mn(mu, data, ind)
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
                1 0 -1;
                0 0 0]; % 18 neighbours + itself

            C.S = spalloc(size(data.w,1), size(data.w,1), size(data.w,1)*19);

            load('s.mat');

            C.S = q;
            return;
            

            for i = 1:data.dim(1)
                for j = 1:data.dim(2)
                    for k = 1:data.dim(3)

                        ind0 = sub2ind(data.dim(1:3), i, j, k);
                        
                        if (data.w(ind0) < eps), continue; end

                        mind = ...
                            ((dind(:,1) + i) > 0) & ...
                            ((dind(:,2) + j) > 0) & ...
                            ((dind(:,3) + k) > 0) & ...
                            ((dind(:,1) + i) <= data.dim(1)) & ...
                            ((dind(:,2) + j) <= data.dim(2)) & ...
                            ((dind(:,3) + k) <= data.dim(3));

                        if (sum(mind) < 4), continue; end

                        c = [i j k] + dind(mind,:);

                        1;

                        ind  = sub2ind(data.dim(1:3), c(:,1), c(:,2), c(:,3));

                        % expected standard deviation
                        s_std = 10;

                        wv = exp(- (data.w(ind) - data.w(ind0)).^2 / (2 * s_std.^2));

                        wv = wv / sum(wv);

                        % [~,ind1] = sort(abs(data.w(ind0) - data.w(ind)), 1, 'ascend');

                        C.S(ind0, ind) = wv;

                    end
                end
            end

            1;


        end

        function [f,b] = do_iter(C, x)

            x_tmp = x;
            
            if (~isempty(C.u))
                x_tmp = x_tmp + C.u; % see Eq. 10 little engine
            end

            C.z = C.image_filter(x_tmp, C.ind);

            [f,b] = C.admm_update(x, C.z);

        end

        function x_flt = image_filter(C, x, ind)

            x_flt = x.new(zeros(size(x.w)));
            x_flt.w = C.S * x.w;
            
        end

    end
end