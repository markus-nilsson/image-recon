classdef data_obj_from_nii_fn < data_obj_image_volume

    properties
        nii_fn;
    end

    methods

        function O = data_obj_from_nii_fn(nii_fn)
            [I,h] = mdm_nii_read(nii_fn);
            O = O@data_obj_image_volume(double(I), h);
            O.nii_fn = nii_fn;
        end
    end
end

