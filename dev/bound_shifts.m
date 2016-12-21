function X = bound_shifts(Y,th)

sz = size(Y);
Mx = spdiags(ones(sz(2),1)*[-1,1],[-1,0],sz(2),sz(2)-1);
My = spdiags(ones(sz(1),1)*[1,-1],[0,1],sz(1)-1,sz(1));

if  max(max(abs(My*Y))) > th || max(max(abs(Y*Mx))) > th
    cvx_begin quiet
        variable X(size(Y))
        minimize norm(X(:)-Y(:),2)
        subject to 
            max(max(abs(My*X))) <= th;
            max(max(abs(X*Mx))) <= th;
    cvx_end
else
    X = Y;
end