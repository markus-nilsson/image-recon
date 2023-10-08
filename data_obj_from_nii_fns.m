classdef data_obj_from_nii_fns < data_obj_image_cell

    properties

        fns = {};

    end

    methods

        function O = data_obj_from_nii_fns(fns)
            assert(iscell(fns), 'input must be cells with file names');

            O = O@data_obj_image_cell(numel(fns));
            O.fns = fns;

            for c = 1:numel(fns)
                O.data_obj{c} = data_obj_from_nii_fn(fns{c});
            end

        end

        function O = trim(O, ir, jr, kr, vr)
            for c = 1:numel(O.data_obj)
                O.data_obj{c}.trim(ir, jr, kr, vr);
            end
        end

    end

end