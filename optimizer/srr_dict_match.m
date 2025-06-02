function [m, ind,w] = imr_dict_match(dict, coeff, c_type)
% function [m, ind] = imr_dict_match(dict, coeff)
%
% 

% xxx: scale free matching? for later
if (nargin < 3), c_type = 3; end

switch (c_type)
    case 1
        TMP = sum( (coeff - dict.W).^2, 2);
        
        [~,ind] = min(TMP);
    
    case 2
        TMP = corr(coeff, dict');

        if (max(TMP) < 0.4)
            1;
        end

        w = exp(- (max(TMP)-TMP).^2 / (2 * 0.01)^2 );
        w = w / (sum(w) + eps);

        [~,ind] = max(w);

        m = w * dict;
        
        m = m / (sqrt(sum(m.^2)) + eps);
        
    case 3
        
        TMP = dict * coeff;

        TMP = TMP ./ sqrt(sum(dict.^2,2));

        [~,ind] = max(TMP);

        m = dict(ind,:);

        % normalize it
        m = m / (sqrt(sum(m.^2)) + eps);
        
        
end


