% Supplementary function for Convergent Cross Mapping Code - Thomas Merkh, tmerkh@ucla.edu - July 25th, 2017
function R = corrcoef(A , B)
  R = [corr(A(:), A(:)) corr(A(:), B(:)); corr(B(:), A(:)) corr(B(:), B(:))];
endfunction