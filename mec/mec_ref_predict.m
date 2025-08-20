function data_ref = mec_ref_predict(data, opts)

if (nargin < 2), opts.present = 1; end

opts = msf_ensure_field(opts, 'lambda', 0.15);
opts = msf_ensure_field(opts, 'c_case', 1);

switch (opts.c_case)

    case 1 % brain tissue
        m = 1e6;

        % Directions
        n = randn(m, 3);
        n = n ./ sqrt(sum(n.^2,2));

        % White matter
        ad = rand(m, 1) * 1.5e-9 + 1.5e-9; % 1.5-3.0 um2/ms
        rd = rand(m, 1) * 0.5e-9 + 0.2e-9; % 0.2-0.7 um2/ms

        dt_wm = tm_1x3_to_1x6(ad, rd, n);

        % Gray matter
        ad = rand(m, 1) * 0.5e-9 + 0.6e-9;
        rd = rand(m, 1) * 0.5e-9 + 0.6e-9;

        dt_gm = tm_1x3_to_1x6(ad, rd, n);

        % CSF
        n = n(1:round(0.1*end), :);
        m = size(n, 1);

        ad = rand(m, 1) * 0.3e-9 + 2.9e-9;
        rd = rand(m, 1) * 0.3e-9 + 2.9e-9;

        dt_csf = tm_1x3_to_1x6(ad, rd, n);

        % signal
        S = cat(2, ...
            exp(-data.xps.bt * dt_wm'), ...
            exp(-data.xps.bt * dt_gm'), ...
            exp(-data.xps.bt * dt_csf'));

        % Some partial volumes and crossing fibers?


    case 2 % Shell grouping (testing needed here)

        S = exp(-data.xps.b * linspace(0.5, 3.0, 200) * 1e-9);

end


c_method = 1;
switch (c_method)

    case 1
        % Build references using approach similar to a Gaussian Process

        K = S * S' / size(S,2);

        K2 = inv(K + opts.lambda * eye(size(K,1))) * K;

        data_ref = data.new( data.w * K2);

        % Rescale
        sc = diag(data.w' * data_ref.w) ./ diag(data.w' * data.w);

        data_ref.w = data_ref.w * diag(1./sc);

        data_ref.w(data_ref.w(:) < 0) = 0;

        1;

    case 2

        % Build references using PCA-like denoising approach (worse)
        if (0) % may be a better target in later runs
            [a,b] = eigs(K, size(K,1));
            K3 = a * diag( 1:size(K,1) < 20 ) * a'; % cut higher eigenvalues

            data_ref = data.new( max(0, data.w * K3));
        end
end