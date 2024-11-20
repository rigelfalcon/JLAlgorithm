% Old method:

p=13;
S1=0;
for k=1:2^p
    S1=S1+k^3;
end
% S1= 1.126174801526784e+15


% New method:
S2=0;
for k=1:p
    S2=S2+2^(p-k)*8^k;
end
% S2=7.330077409280000e+11


% The improvement is rather big as you see.

Lasha