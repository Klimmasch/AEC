function y  = laprnd(m, n, mu, sigma, boarders)
    %LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
    %   with mean mu and standard deviation sigma. 
    %   mu      : mean
    %   sigma   : standard deviation
    %   [m, n]  : the dimension of y.
    %   boarders: min and max values for y.
    %   Default mu = 0, sigma = 1. 
    %   For more information, refer to
    %   http://en.wikipedia.org./wiki/Laplace_distribution

    %   Author  : Elvis Chen (bee33@sjtu.edu.cn)
    %   Date    : 01/19/07

    %Check inputs
    if nargin < 2
        error('At least two inputs are required');
    end

    if nargin == 2
        mu = 0; sigma = 1;
        boarders = [];
    end

    if nargin == 3
        sigma = 1;
        boarders = [];
    end

    if ~isempty(boarders)
        y = ones(m, n) * Inf;
        iters = 1;
        while (min(min(y))<boarders(1)) || (max(max(y))>boarders(2))
            % Generate Laplacian noise
            u = rand(m, n)-0.5;
            b = sigma / sqrt(2);
            y = mu - b * sign(u).* log(1- 2* abs(u));
            iters = iters + 1;
            if iters > 1000
                error('You might want to change the boarder conditions, bro.')
            end
        end
    else
        u = rand(m, n)-0.5;
        b = sigma / sqrt(2);
        y = mu - b * sign(u).* log(1- 2* abs(u));
    end
end