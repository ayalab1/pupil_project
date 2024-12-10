% smoothed = SMOOTHVECT(origsignal, filter)
%            convolves origsignal with filter, removes length(filter) - 1
%            points to produce a signal with the same length as origsignal, and
%            returns the result
function [s] = smoothvect(os, f)

sm = conv(f, os);
filterlen = length(f);
% the convolution returns a vector whose length is the sum of
% the length of the two arguments - 1
startp = round(filterlen / 2);
endp = startp + length(os) - 1;
s = sm(startp:endp);
