function X=fast_hankel(x)
t=0:size(x,1)-1;
[tx, ty] = ndgrid(t);

xf = griddedInterpolant(t, x, 'nearest', 'none');

X=xf(tx + ty);
X(isnan(X)) = 0;

end
