classdef mec_params_slice < mec_params_base

    methods

        function o = mec_params_slice(data)
            o = o@mec_params_base(data);

            o.opt_opts.FunctionTolerance = 1e-4;

        end

        function x = get_x(o, c_vol)
            x = o.x(:,:,c_vol);
            x = x(:);
        end

        function p = get_p(o, c_vol) % scale later
            p = o.get_x(c_vol);
            p = reshape(p, 3, numel(p) / 3);
            
            p(1,:) = p(1,:);
            p(2,:) = p(2,:) / 100;
            p(3,:) = p(3,:) / 100;
        end

        function o = update(o, x, c_vol)
            o.x(:,:,c_vol) = reshape(x, size(o.x, [1 2]));
            o.enabled = true;
        end

        function o = reset(o)
            o.x = zeros(3, o.data.dim(3), o.data.n_vol);
            o.enabled = false;
        end


        function c = distort(o, c, c_vol)

            if (~o.enabled), return; end

            p = o.get_p(c_vol);

            % Displacement, scale, and shear
            displ_y = 0 + p(1, c.z_ind)';
            scale_y = 1 + p(2, c.z_ind)';
            shear_y = 0 + p(3, c.z_ind)';

            c.y = scale_y .* c.y + shear_y .* c.x + displ_y;
            
        end

        function plot(o)

            p = o.get_p(1:o.data.n_vol);

            msf_clf;
            yscl = [1 0.05 0.05];
            
            for c_param = 1:3

                subplot(3,3,c_param);
                plot(squeeze(p(c_param,:,:)));
                ylim([-1 1] * yscl(c_param));
                switch (c_param)
                    case 1
                        ylabel('Displacement y [vox]');
                    case 2
                        ylabel('Scale y');
                    case 3
                        ylabel('Shear y');
                end
                xlabel('Slice');
                xlim([0 o.data.dim(3)+1]);

                subplot(3,3,c_param + 3);
                plot(squeeze(mean(abs(p(c_param,:,:)), 3)));
                ylim([-1 1] * yscl(c_param));
                switch (c_param)
                    case 1
                        ylabel('mean(abs(dy)) [vox]');
                    case 2
                        ylabel('Scale y');
                    case 3
                        ylabel('Shear y');
                end
                xlabel('Slice');
                xlim([0 o.data.dim(3)+1]);  


                % regression
                try
                    X = [ones(o.data.n_vol, 1)';
                        o.u(:,1)';
                        o.u(:,2)';
                        o.u(:,3)';
                        o.v(:,1)';
                        o.v(:,2)';
                        o.v(:,3)';
                        ];

                    for c_slice = 1:size(p,2)
                        tmp = squeeze(p(c_param, c_slice, :));
                        beta(:, c_slice) = tmp' * X' * inv(X * X');

                        tmp_hat = beta(:,c_slice)' * X;

                        SS_res = sum((tmp(:) - tmp_hat(:)).^2);              % residual sum of squares
                        SS_tot = sum((tmp - mean(tmp)).^2);            % total sum of squares
                        R2(c_slice) = 1 - SS_res/SS_tot;
                    end

                    if (1)
                        subplot(3, 3, 6 + c_param)
                        plot(R2);
                        ylim([0 1]);
                        xlim([0 o.data.dim(3)+1]);

                        ylabel('R2 of g(n, n-1)');
                    else

                        x = 1:size(beta, 2);
                        for c = 1:3
                            subplot(3, 9, 18 + c + (c_param-1) * 3);
                            ind = 1:1:(numel(x));
                            plot(x(ind), beta(c+1, ind), x(ind), beta(c+4, ind), 'linewidth', 1);


                            axis off
                        end
                    end
                catch me
                end

            end

            pause(0.05);

        end

    end
end
