function g = orientedGaussian(fSize, sigma1, sigma2)
    ind = -floor(fSize/2) : floor(fSize/2);
    g1 = exp(-(ind.^2) / (2*sigma1^2));
    g2 = exp(-(ind.^2) / (2*sigma2^2));
    g = g1' * g2;
    g = g / sum(g(:));
end