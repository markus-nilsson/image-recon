classdef do < handle

    properties
        transpose = 0;
        c_type; % 1 - weights; 2 - cell
    end

    methods

        % a do (data object) holds either the acquired data volume(s) or 
        % the estimated image / model coefficients
        % 
        % we can describe it as w_il, where i refers to voxel indices and 
        % l to the index of the experimental setting (e.g. b-value/direction
        % or more abstract such as phase encoding direction
        %
        % the data object has two embodiments, one where the weights are 
        % stored as a matrix (w_il), and one where the second dimension
        % is stored as cells (w_i){l}. in the latter case, the number of 
        % voxels (i = 1..n_vox) can depend on l, such that n_vox(l)

        function O = do()
        end

        function O = ctranspose(O)
            O.transpose = ~O.transpose;
        end

    end

    methods (Abstract)

%         O_new = new(O, w)
%         r = norm(O)
%         s = sum(O)
%         n = numel(O)
% 
%         O_new = add_noise_gaussian(O, noise_std)

    end

end

