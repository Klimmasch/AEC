function [h, array] = display_network(A, opt_normalize, opt_graycolor, cols, opt_colmajor, opt_zeromean)
% This function visualizes filters in matrix A. Each column of A is a
% filter. We will reshape each column into a square image and visualizes
% on each cell of the visualization panel.
% All other parameters are optional, usually you do not need to worry
% about it.
% opt_normalize: whether we need to normalize the filter so that all of
% them can have similar contrast. Default value is true.
% opt_graycolor: whether we use gray as the heat map. Default is true.
% cols: how many columns are there in the display. Default value is the
% squareroot of the number of columns in A.
% opt_colmajor: you can switch convention to row major for A. In that
% case, each row of A is a filter. Default value is false.
% opt_zeromean: rescale input patches to zero mean. Default value is true.
warning off all

if ~exist('opt_normalize', 'var') || isempty(opt_normalize)
    opt_normaliz e= true;
end

if ~exist('opt_graycolor', 'var') || isempty(opt_graycolor)
    opt_graycolor = true;
end

if ~exist('opt_colmajor', 'var') || isempty(opt_colmajor)
    opt_colmajor = false;
end

if ~exist('opt_zeromean', 'var') || isempty(opt_zeromean)
    opt_zeromean = true;
end

% rescale
if opt_zeromean
    A = A - mean(A(:));
end

if opt_graycolor, colormap(gray); end

% compute rows, cols
[L M] = size(A);
sz = sqrt(L);
buf = 1;
if ~exist('cols', 'var')
    if floor(sqrt(M))^2 ~= M
        n = ceil(sqrt(M));
        while ((mod(M, n) ~= 0) && (n < 1.2 * sqrt(M)))
            n = n + 1;
        end
        m = ceil(M/n);
    else
        n = sqrt(M);
        m = n;
    end
else
    n = cols;
    m = ceil(M/n);
end

array = -ones(buf+m*(sz+buf), buf+n*(sz+buf));

if ~opt_graycolor
    array = 0.1.* array;
end


if (~opt_colmajor)
    k = 1;
    for (i = 1:m)
        for (j = 1:n)
            if (k > M)
                continue;
            end
            clim = max(abs(A(:, k)));
            if (opt_normalize)
                array(buf+(i-1)*(sz+buf)+(1:sz), buf+(j-1)*(sz+buf)+(1:sz)) = reshape(A(:, k), sz, sz)/clim;
            else
                array(buf+(i-1)*(sz+buf)+(1:sz), buf+(j-1)*(sz+buf)+(1:sz)) = reshape(A(:, k), sz, sz)/max(abs(A(:)));
            end
            k = k + 1;
        end
    end
else
    k = 1;
    for (j = 1:n)
        for (i = 1:m)
            if (k > M)
                continue;
            end
            clim = max(abs(A(:, k)));
            if (opt_normalize)
                array(buf+(i-1)*(sz+buf)+(1:sz), buf+(j-1)*(sz+buf)+(1:sz)) = reshape(A(:, k), sz, sz)/clim;
            else
                array(buf+(i-1)*(sz+buf)+(1:sz), buf+(j-1)*(sz+buf)+(1:sz)) = reshape(A(:, k), sz, sz);
            end
            k = k + 1;
        end
    end
end

%rescales ouputs for plot
if (opt_zeromean)
    lims = [-1 1];
end

if (opt_graycolor)
    h = imagesc(array,'EraseMode','none', lims);
else
    h = imagesc(array,'EraseMode','none', lims);
end
axis image off

drawnow;

warning on all

%add central lines
% hold on
% line([ceil(size(array, 2)/2) ceil(size(array, 2)/2)], [1 ceil(size(array, 1))], ...
%     'Color','r','LineWidth', 3) ; %vertical
% line([1 ceil(size(array, 2))], [ceil(size(array, 1)/2) ceil(size(array, 1)/2)] , ...
%     'Color','r','LineWidth', 3) ; %hor
% hold off
