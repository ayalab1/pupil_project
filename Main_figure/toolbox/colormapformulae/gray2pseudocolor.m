% It converts gray image to pseudocolor image
%
% Inputs:
%  x: gray image whose range is [0,1]
%  type: colormap type
%
% Output:
%  col: pseudo color image [0,1]
%
% typical colormap names of pnm3d:
%  [ 7, 5,15]:traditional,pm3d,black-blue-red-yellow
%  [ 3,11, 6]:green-red-violet
%  [23,28, 3]:ocean,green-blue-white
%  [21,22,23]:hot,black-red-yellow-white
%  [30,31,32]:color printable on gray,black-blue-violet-yellow-white
%  [33,13,10]:rainbow,blue-green-yellow-red
%  [34,35,36]:AFM hot,black-red-yellow-white
%
% functions:
%   0: 0
%   1: 0.5
%   2: 1
%   3: x
%   4: x^2
%   5: x^3
%   6: x^4
%   7: sqrt(x)
%   8: sqrt(sqrt(x))
%   9: sin(90x)
%  10: cos(90x)
%  11: |x-0.5|
%  12: (2x-1)^2
%  13: sin(180x)
%  14: |cos(180x)|
%  15: sin(360x)
%  16: cos(360x)
%  17: |sin(360x)|
%  18: |cos(360x)|
%  19: |sin(720x)|
%  20: |cos(720x)|
%  21: 3x
%  22: 3x-1
%  23: 3x-2
%  24: |3x-1|
%  25: |3x-2|
%  26: (3x-1)/2
%  27: (3x-2)/2
%  28: |(3x-1)/2|
%  29: |(3x-2)/2|
%  30: x/0.32-0.78125 
%  31: 2*x-0.84
%  32: 4x;1;-2x+1.84;x/0.08-11.5
%  33: |2*x - 0.5|
%  34: 2*x
%  35: 2*x - 0.5
%  36: 2*x - 1
% negative numbers mean inverted (1-output)
%
% version:201800402


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colormap formulae:                                       %
%                                                          %
% Copyright (C) 2018 Masayuki Tanaka. All rights reserved. %
%                    mtanaka@sc.e.titech.ac.jp             %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function col = gray2pseudocolor( x, type )

if( ~exist('type','var') || isempty(type) )
 type = [7,5,15];
end

if( isa(type, 'char') )
 type = type2ids( type );
end

col = zeros( size(x,1), size(x,2), 3 );

for i=1:3
 col(:,:,i) = formulae( x, type(i) );
end

col(col<0) = 0;
col(col>1) = 1;

end