classdef mec_params_volume_rigid < mec_params_base
    
    methods

        function o = mec_params_volume_rigid(data)
            o = o@mec_params_base(data);
        end

        function o = reset(o)
            for c_vol = 1:o.data.xps.n
                o = o.update(zeros(1,6), c_vol);
            end
            o.enabled = false;
        end

        function o = update(o, x, c_vol)
            o.x(:,c_vol) = x;
            o.enabled = true;
        end

        function x = get_x(o, c_vol)
            x = o.x(:,c_vol);
        end        

        function p = get_p(o, c_vol)
            x = o.get_x(c_vol);
            p = x / 100;
        end

        function d = distort(o, c, c_vol)

            if (isempty(c)) || (~o.enabled)
                d = c;
                return; 
            end
            
            % Apply volume-wise corrections
            p = o.get_p(c_vol);

            % Translations...
            c.x = c.x + p(1);
            c.y = c.y + p(2);
            c.z = c.z + p(3);

            d = c; % Pull in other parameters            

            % ...and rotations
            R = mec_params_volume.euler_xyz(p(4), p(5), p(6));
            
            d.x = R(1,1) * c.x + R(1,2) * c.y + R(1,3) * c.z;
            d.y = R(2,1) * c.x + R(2,2) * c.y + R(2,3) * c.z;
            d.z = R(3,1) * c.x + R(3,2) * c.y + R(3,3) * c.z;

        end
        

        function plot(o)

            p = o.get_p(1:o.data.n_vol)';

            p(:,4:6) = p(:,4:6) / pi * 180;

            t_str = {'dx', 'dy', 'dz', 'rx', 'ry', 'rz'};
            
            msf_clf;

            yl = [1 2 1  3 3 3];
            for c = 1:6

                subplot(2,4,c + (c > 3));
                plot(p(:,c));
                title(t_str{c});
                ylim([-1 1] * yl(c));
                
                if (c <= 3)
                    ylabel('Voxels [1]');
                else
                    ylabel('Degrees [1]');
                end

            end

            pause(0.05);
            
        end

        function o = show(o)
            error('not updated');
            J = [];
            for c = 1:6

                switch (c)

                    case 1
                        I0 = data_mov.imreshape(c_vol);
                        I = I0;
                    case 2
                        ttmp = mec_points_all(data_mov).apply(params);
                        ID = ttmp.imreshape(c_vol);
                        I = ID;
                    case 3
                        IR = data_ref.imreshape(c_vol);
                        R0 = I0 - IR;
                        R1 = I  - IR;
                        I = IR;
                    case 4
                        I = R1; % registered
                    case 5
                        I = R0; % original
                    case 6
                        I = (R0 - R1);
                        I = I / (max(abs(I(:))) / max(abs(R0(:))) + eps);
                end

                k = round([0.15 0.5 0.85] * data_mov.dim(3));
                J = cat(2, cat(1, I(:,:,k(1)), I(:,:,k(2)), I(:,:,k(3))), J);

                if (c == 3)
                    subplot(1,2,1);
                    msf_imagesc(J);
                    J = [];
                elseif (c == 6)
                    subplot(1,2,2);
                    msf_imagesc(J);
                    caxis([-1 1] * 0.6 * (max(0.1, max(abs(J(:))))));
                    1;
                end

            end

            pause(0.05);
        end
    end
        

    methods (Static)

        function [R, Rinv] = euler_xyz(alpha, beta, gamma)
            % alpha = rotation about x
            % beta  = rotation about y
            % gamma = rotation about z
            %
            % Active extrinsic XYZ (roll, pitch, yaw)
            % Use radians.

            Rx = [1 0 0;
                0 cos(alpha) -sin(alpha);
                0 sin(alpha) cos(alpha)];

            Ry = [cos(beta) 0 sin(beta);
                0 1 0;
                -sin(beta) 0 cos(beta)];

            Rz = [cos(gamma) -sin(gamma) 0;
                sin(gamma)  cos(gamma) 0;
                0           0          1];

            R = Rz * Ry * Rx;
            Rinv = R.'; % exact inverse
        end

    end
end
