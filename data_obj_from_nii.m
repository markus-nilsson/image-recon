classdef data_obj_from_nii < data_obj_cell

    properties
       

    end

    properties (Access = protected)

    end

    properties (Access = private)

    end

    methods

        function O = data_obj_from_nii_fn(nii_fn)
            % assume input to be cell of nii filenames
            

             for c = 1:numel(fns)
                 O.data{c} = data_obj_from_niis(nii_fns{c});
             end

        end


    end
end

