classdef do_w_from_nii < do_w_image_volume

    properties
        nii_fn;
        xps = [];
    end

    methods

        function O = do_w_from_nii(nii_fn)
            [I,h] = mdm_nii_read(nii_fn);
            O = O@do_w_image_volume(double(I), h);
            O.nii_fn = nii_fn;


            try 
                O.xps = mdm_xps_from_nii_fn(nii_fn, 1);
            catch
            end
        end
    end
end

