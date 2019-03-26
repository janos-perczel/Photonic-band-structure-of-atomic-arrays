function [ zComplexOut ] = PosSqrt( zComplexIn )
%PosSqrt (Positive Square Root) takes a complex number, z, as its input and returns
%its complex square root, where both the real and imaginary part are
%greater than (or equal to) zero. This is consistent with the so-called Sommerfeld
%integration path IN THE FIRST QUADRANT that ensures that the radiation condition is satisfied.

%NB: This choice of square root is valid only in the region on the top Riemann
%sheet of z that corresponds to the first quadrant of the sqrt(z) plane.

zSqrt=sqrt(zComplexIn);

zSqrtReal=real(zSqrt);
zSqrtReal=abs(zSqrtReal);

zSqrtImag=imag(zSqrt);
zSqrtImag=abs(zSqrtImag);

zComplexOut=zSqrtReal+1i*zSqrtImag;

end

