%%% Generation of a set of Gaussian filters
% @param fSize              total filter size
% @param sigma1             variance in vertical direction
% @param sigma2             variance in horizontal direction
%%%
function [g1, g2] = orientedGaussianVectors(fSize, sigma1, sigma2)
    ind = -floor(fSize/2) : floor(fSize/2);
    g1 = exp(-(ind.^2) / (2*sigma1^2));
    g1 = g1 / sum(g1);
    g2 = exp(-(ind.^2) / (2*sigma2^2));
    g2 = g2 / sum(g2);
    % g = g1' * g2; % 2D-filter, if necessary
    % g = g / sum(g(:));
end
