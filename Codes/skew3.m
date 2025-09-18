%%HELPER FUNCTION (MEKF)
function S = skew3(v)
S = [  0   -v(3)  v(2);
      v(3)   0   -v(1);
     -v(2) v(1)    0 ];
end
