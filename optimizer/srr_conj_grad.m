function [x,cgr] = srr_conj_grad(E,b,x,opts)
% A conjugate gradient solver that allows to plot the result in each
% iteration
% Solves: E*x=b
%
% x = conjugate_gradient(E,b)
% x = conjugate_gradient(E,b,tol)
% x = conjugate_gradient(E,b,tol,maxit)
% x = conjugate_gradient(E,b,tol,maxit,z)
% x = conjugate_gradient(E,b,tol,maxit,z,verbose)
% x = conjugate_gradient(E,b,tol,maxit,z,verbose,newfigure)
% [x,x_iter] = conjugate_gradient(_____)
%
%
% Input:
%   E         =  Hermitian and positive definite matrix or operator that
%                immitates such a matrix and implements a mtimes funciton.
%                Alternatively, A can function handle that imitates such a
%                matrix.
%   b         =  Right hand side of the equation A*x=b
%   tol       =  Tolerance of the method (default: 1e-6)
%   maxit     =  maximum number of iterations (default: 20)
%   x         =  initial guess (default: 0)
%   verbose   =  0 for no output, 1 for plotting the images and in each
%                iteration and 2 for printing the residal of the for each
%                iteration.
%   newfigure =  0 for using the same figure for each call of the funciton
%                or 1 for opening a new figure.
%
%
% Output:
%   x      = Result
%   x_iter = Result after each iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col = @(x) x(:);
% norm = @(x) sum(col(x).^2);

if (1)
    Ed = E(x);
    d = b - Ed;
    r = d;
    
    normrr0 = norm(b);%.norm();
    normrr  = norm(r);%r.norm();
    
    for n = 1:opts.n_iter_cg
        
        Ed = E(d);
        tmp = conj(d) .* Ed;

        alpha = normrr/real(sum(col(tmp))+eps);
        x = x + d * alpha;

        r = r - Ed * alpha;
        normrr2 = norm(r) + eps;%norm();
        beta = normrr2/normrr;
        normrr = normrr2;
        d = r + d * beta;

        if (numel(opts.cg_display_ind) > 0) && (mod(n-1, 1) == 0)
            msf_clf;
            tmp = x.imreshape();

            switch (numel(opts.cg_display_ind))
                case 1
                    nr = 1;
                    nc = 1;
                case 2
                    nr = 1;
                    nc = 2;
                case 3
                    nr = 1;
                    nc = 3;
                case 4
                    nr = 2;
                    nc = 2;
                case 5
                    nr = 1;
                    nc = 5;
                case 6
                    nr = 2;
                    nc = 3;
                otherwise
                    nr = 2;
                    nc = 3;
                    opts.cg_display_ind = opts.cg_display_ind(1:6);
            end

            for c = 1:numel(opts.cg_display_ind)
                subplot(nr,nc,c);
                msf_imagesc(real(tmp(:,:,:,opts.cg_display_ind(c))));
                colorbar;
                title(num2str(sqrt(normrr/normrr0)));
            end
            pause(0.005);

        end

        cgr(n) = sqrt(normrr/normrr0);

        if (cgr(n) < opts.cg_tol), break; end
        
    end

end
