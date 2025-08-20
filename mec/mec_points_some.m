classdef mec_points_some < mec_points

    methods

        function o = mec_points_some(data, mask, frac)
            o = o@mec_points(data);

            % select random points in mask
            ind = find(mask(:));
            k = round( frac * sum(mask(:)) );

            o.ind = ind(randperm(numel(ind), k));

            o.coords = o.ind2coords(o.ind);

        end

        function show_mask(o, mask)
            % 
            % % Plot
            % figure(my_fig('mec_points_some'));
            % 
            % k_list = round(linspace(0.1, 0.9, 9) * data.dim(3));
            % K = []; J = [];
            % for c_k = 1:numel(k_list)
            % 
            %     J = cat(1, J, mask(:,:,k_list(c_k)));
            % 
            %     if (mod(c_k - 1, 3) == 2)
            %         K = cat(2, K, J);
            %         J = [];
            %     end
            % 
            % 
            % end
            % 
            % subplot(1,2,1);
            % msf_imagesc(K);
            % caxis([0 1]);
            % colorbar;


        end

        function ind = get_ind(o)
            ind = o.ind;
        end

        function y = apply(o, data, params, c_vol)

            if (nargin < 3), params = []; end
            if (nargin < 4), c_vol = []; end


            if (isempty(c_vol))

                for c_vol = data.n_vol:-1:1
                    y(:,c_vol) = o.apply_one_vol(data, c_vol, params);
                end

                y = y(:);

            else

                y = o.apply_one_vol(data, c_vol, params);

            end

        end

    end

end