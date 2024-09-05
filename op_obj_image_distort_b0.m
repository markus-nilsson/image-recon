classdef op_obj_image_distort_b0 < op_obj_image_distort

    methods

        function O = op_obj_image_distort_b0(B0, I, pe, n_k)

            % for now, assume B0 and h have the same FOV, and that the
            % phase encoding is along the second dimension of B0
            % thus, assume here we shift in the y direction

            t = double(B0.imreshape()) * pe * B0.pixdim(2);

            z = zeros(numel(t), 1);
            t = cat(2, z, t(:), z);

            d = I.new(t);

            O = O@op_obj_image_distort(d, n_k);

        end

    end
end