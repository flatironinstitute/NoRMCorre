function y = han(x,frac)
    % apply a hamming window to signal x of length x*(1+frac).
    % equivalent to zero-padding, windowing and then truncating.
    
    if nargin < 2 || isempty(frac)
        frac = 0.5;
    end
    
    if isinf(frac); 
        y = x;
    else
        sx = size(x);
        h1 = hamming(sx(1) + 2*round(sx(1)*frac));
        h2 = hamming(sx(2) + 2*round(sx(2)*frac));
        %y = x.*(h1(round(sx(1)*frac)+(1:sx(1)))*h2(round(sx(2)*frac)+(1:sx(2)))');
        y = bsxfun(@times, bsxfun(@times,h2(round(sx(2)*frac)+(1:sx(2)))',x), h1(round(sx(1)*frac)+(1:sx(1))));
        if length(sx) == 3
            temp = hamming(sx(3) + 2*round(sx(3)*frac));
            h3(1,1,1:sx(3)) = temp(round(sx(3)*frac)+(1:sx(3)));
            y = bsxfun(@times, y, h3);
        end            
    end