% It return mx3 matrix which represents colormap.
%
% Inputs:
%  type: name of colormap or
%        three dimensional function id
%  m: number of data (default: 256)
%
% Output:
%  c: mx3 sized colormap matrix
%
% Example:
%  % surf();
%  % colormap(colormapformulae([7,5,15]));
%
% buildin colormap names:
%  parula, jet, hsv, hot, cool, spring, summer, autumn, winter, gray, bone, copper, pink
% % doc colormap % for more information
%
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
function c = colormapformulae( type, m )

if( ~exist('type','var') || isempty(type) )
 type = [7,5,15];
end

if( ~exist('m','var') || isempty(m) )
 m =256;
end

if( strcmpi(type, 'parula') )
 c = parula(m);
elseif( strcmpi(type, 'jet') )
 c = jet(m);
elseif( strcmpi(type, 'hsv') )
 c = hsv(m);
elseif( strcmpi(type, 'cool') )
 c = cool(m);
elseif( strcmpi(type, 'spring') )
 c = spring(m);
elseif( strcmpi(type, 'autumn') )
 c = autumn(m);
elseif( strcmpi(type, 'winter') )
 c = winter(m);
elseif( strcmpi(type, 'bone') )
 c = bone(m);
elseif( strcmpi(type, 'copper') )
 c = copper(m);
elseif( strcmpi(type, 'pink') )
 c = pink(m);
 
else
 if( isa(type, 'char') )
  type = type2ids( type );
 end

 x = [0:1.0/(m-1):1];
 c = zeros(numel(x),3);
 for i=1:3
  c(:,i) = formulae(x,type(i));
 end
 c(c<0)=0;
 c(c>1)=1;
end

end
