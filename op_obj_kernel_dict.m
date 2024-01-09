classdef op_obj_kernel_dict < op_obj_kernel

    properties
        S_k;
        CEk;
        mp_list;

        Mw;

    end
    
    methods
        
        function O = op_obj_kernel_dict(S_k, mp_list, n_comp, n_vox)

            if (numel(mp_list) ~= (ndims(S_k)-1))
                error('expecting a %i cells for mp_list', ndims(S_k)-1);
            end

            S_k = reshape(S_k, [size(S_k,1), numel(S_k)/size(S_k,1)]);


            % This works better, although I think it would work worse
            if (1) % high values in COEFF (which we will estimate)


                [COEFF, SCORE] = pca(S_k, ...
                    'centered', 'off', ...
                    'numcomponents', n_comp); % auto-select this

            else % low values in COEFF, high in SCORE

                [COEFF2, SCORE2] = pca(S_k', ...
                    'centered', 'off', ...
                    'numcomponents', n_comp); % auto-select this
                
                SCORE = COEFF2;
                COEFF = SCORE2;
                clear SCORE2 COEFF2;

            end

            % Reduce number of components if possible
            for c = 1:n_comp
                R(c) = mean(sqrt(mean( (S_k - SCORE(:,1:c) * COEFF(:,1:c)').^2)));
            end

            n_comp = sum(R > 1e-6);

            SCORE = SCORE(:, 1:n_comp);
            COEFF = COEFF(:, 1:n_comp);

            % normalize the coefficients to unit length
            COEFF = COEFF ./ sqrt(sum(COEFF.^2,2));

            O = O@op_obj(SCORE);


            O.CEk = COEFF;
            O.S_k = S_k;
            O.mp_list = mp_list;

            O.Mw = sqrt(diag(O.M' * O.M));


        end

        function x = init_x(O, c_mode)

            if (nargin < 2), c_mode = 1; end

            switch (c_mode)
                case 1
                    x = O.init_zero();
                case 2
                    x = repmat(mean(O.CEk,1)', [1 O.n_vox]);
                    x  = x .* mean(data ./ (O * x),1);
            end
        end

        function x = init_zero(O)
            x = zeros(O.n_mp, O.n_vox);
        end

        function x = init_mp_with_db(O, s)

            tmp = O.DBmerge' * s;

            tmp = O.S_k' * tmp;

            tmp = tmp ./ sqrt(sum(O.S_k.^2,1))';

            [~,ind] = max( tmp , [], 1 );
            
            x = O.CEk(ind,:)';
            
            x  = x .* mean(data ./ (O * x),1);

        end        
        
        function D = dict(O)
            D = O.CEk;
        end


        function d = map_to_dict(O, x)
            
            Q = O.dict(); 

            % assumes a normalized dictionary
            PP = Q * diag(sqrt(O.Mw)) * x;
            [~,ind] = max(PP,[],1);

            d = Q(ind,:)';

        end

        function m = x_to_m(O, x)

            m = zeros(numel(O.mp_list), size(x,2));
            Q = O.dict();

            sz = zeros(size(O.mp_list));
            for i = 1:numel(O.mp_list)
                sz(i) = numel(O.mp_list{i});
            end

            for i = 1:numel(O.mp_list)
                tmp_sz = sz;
                tmp_sz(i) = 1;
                order = 1:(2+numel(O.mp_list));
                order([2 2+i]) = order([2+i 2]);
                A{i} = squeeze(permute(repmat(O.mp_list{i}, [1 1 tmp_sz]), order));
            end



            for c = 1:size(x,2)
                disp(c/size(x,2));

                tmp = x(:,c);

                tmp = tmp / sqrt(sum(tmp.^2));

%                 w = sum( (Q - tmp').^2, 2);

                % account for effect on signal
                w = (Q - tmp').^2 * (O.Mw).^2;

                % w2 = sum( ( (Q-tmp') * O.M' ).^2, 2);

                [w,ind] = sort(w);

                nw = 2^numel(O.mp_list);

                WX = diag(O.Mw);
%                 WX = O.M;

                w2 = lsqnonneg( (Q(ind(1:nw),:) * WX')', (WX * tmp));
%                 w2 = w2 / sum(w2);

                for i = 1:numel(O.mp_list)
                    tmp2 = A{i}(ind(1:nw));
                    m(i, c) = sum(tmp2 .* w2);
                    if (i == numel(O.mp_list))
                        1;
                    end
                end

            end

        end

        function d = mp_to_dict(O, mp)

%             if (size(mp,2) > 1)
%                 error('not implemented');
%             end
% 
%             for c = 1:size(mp, 1)
% 
%                 sub(c) = find(mp(c) > O.mp_list{c}, 1, 'last');
% 
%                 1;
% 
%                 dcm(c) = (mp(c) - O.mp_list{c}(sub(c))) / ...
%                     (O.mp_list{c}(sub(c)+1) - O.mp_list{c}(sub(c)));
%                 
%             end
% 
%             sz = cell2mat(cellfun(@(x) numel(x), O.mp_list, 'UniformOutput', false));
% 
%             for c = 1:size(mp,1)
% 
%                 ind = sum(sz.^( 0:1:5) .* (sub-1) ) + 1;
% 
%                 d = O.
% 
%                 d(c) = dcm * p2 * (1-dcm) * p1;
% 
%             end
%             
% 
% 
% 
%             1;
% %             ind2sub(size())

        end


        function simplify_dict()

            %             % Simplify the dictionary
            %             n = 3;
            %             idx_i = kmeans(O.CEk, n);
            %             CEkr = zeros(n*n*n, size(O.CEk,2));
            %             c = 1;
            %             for i = 1:n
            %                 CEj = O.CEk(idx_i == i, :);
            %                 idx_j = kmeans(CEj, n);
            %                 for j = 1:n
            %                     CEk = CEj(idx_j == j, :);
            %                     idx_k = kmeans(CEk, n);
            %                     for k = 1:n
            %                         CEkr(c,:) = mean(O.CEk(idx_k == k, :), 1);
            %                         c = c + 1;
            %                     end
            %                 end
            %             end
            %
            %             O.CEk_full = O.CEk;
            %             O.CEk = CEkr;

        end        
        

    end
end
