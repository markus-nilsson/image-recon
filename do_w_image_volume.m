classdef do_w_image_volume < do_w_image_vector
    
    methods

        function O = do_w_image_volume(I, h)

            I = reshape(I, [prod(size(I,1,2,3)) size(I,4)]);
            O = O@do_w_image_vector(I, h);

        end



    end

    methods (Static)

        % not ideal place for this
        function data_obj = new_hr_from_lr(data_lr)

            h_hr = data_lr.h;
            tmp = single(...
                h_hr.dim(4))*h_hr.pixdim(4)/(h_hr.pixdim(2));
            h_hr.dim(4) = int16(tmp);
            h_hr.pixdim(4) = data_lr.h.pixdim(2);

            % also adjust srow's
            h_hr.srow_z(3) = h_hr.pixdim(4);

            X = [h_hr.srow_x'; h_hr.srow_y'; h_hr.srow_z';];
            X = X(1:3, 1:3);
            X = X - diag(diag(X));

            if (sum(X(:).^2) > 1e-5)
                error('need transverse slices for now');
            end

            I = zeros(double(h_hr.dim(2:5)'));

            data_obj = do_w_image_volume(I,h_hr);

        end
    end
end

