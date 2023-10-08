classdef data_obj_image_vector < data_obj

    properties


    end

    properties (GetAccess = public, SetAccess = protected)
        h;
    end

    methods

        function O = data_obj_image_vector(w, h)
            O = O@data_obj(w);
            O.h = h;            
        end

        function O_new = new(O, w)
            assert(prod(O.h.dim(2:4)) == size(w, 1), 'size error');
            assert(prod(O.h.dim(5)) == size(w, 2), 'size error');
            O_new = data_obj_image_vector(w, O.h);
        end

        function O = add(O, w, k)
            if (my_isa(w, 'data_obj')), f = @(x) x.w; else f = @(x) x; end
            O.w(:,k) = O.w(:,k) + f(w);
        end

        function I = imreshape(O)
            I = reshape(O.w, O.h.dim(2:5)');            
        end

        function O_new = zeros(O)
            O_new = data_obj_image_vector(zeros(size(O.w)), O.h);
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

            I = zeros(double(h_hr.dim(2:5)'));

            data_obj = data_obj_image_volume(I,h_hr);

        end
    end
end

