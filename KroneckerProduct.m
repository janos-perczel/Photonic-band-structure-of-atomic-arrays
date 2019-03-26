function tp=KroneckerProduct(varargin)
%takes variable number of inputs and creates their tensor product, from
%left to right. Returns tp as final result.
tp=varargin{1};
for k=1:nargin-1
    tp=kron(tp,varargin{k+1});
end