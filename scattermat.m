function scattermat(pts, marker)
    if (nargin == 1)
        marker = 'b';
    end
    sz = size(pts);
    if(sz(1) == 2)
        scatter(pts(1,:),pts(2,:),marker);
    elseif (sz(2) == 2)
        scatter(pts(:,1), pts(:,2),marker);
    end
end