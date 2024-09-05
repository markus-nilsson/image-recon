classdef op_obj_image_distort_b0h < op_obj_image
    % use high resolution to avoid interpolation artefacts

    properties
        U;
        P;
    end

    methods

        function O = op_obj_image_distort_b0h(B0, I, pe, n_k)

            O = O@op_obj_image([], B0.h, I.h, n_k);

            % Define a HR operator in the y-direction (phase encode)
            h_hr = B0.h;
            sc = 3;
            h_hr.sform_code = 0;
            h_hr.pixdim(3) = h_hr.pixdim(3) / sc;
            h_hr.dim(3) = h_hr.dim(3) * sc;

            O.U = op_obj_image_h2l(B0.h, h_hr, 3);

            B0h = O.U'*B0;

            pe = pe * sc^2;

            O.P = op_obj_image_distort_b0(B0h, B0h, pe, n_k);

            O.aM = O.U.S * O.P.S * O.U.S' / sc;
            O.aMT = O.aM';

        end
    end
end