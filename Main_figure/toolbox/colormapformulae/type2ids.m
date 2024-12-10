% It converts string to ids
%
% Input:
%  type: name
%
% Output:
%  ids: three dimensional id
%
% typical colormap names of pnm3d:
%  [ 7, 5,15]:traditional,pm3d,black-blue-red-yellow
%  [ 3,11, 6]:green-red-violet
%  [23,28, 3]:ocean,green-blue-white
%  [21,22,23]:hot,black-red-yellow-white
%  [30,31,32]:color printable on gray,black-blue-violet-yellow-white
%  [33,13,10]:rainbow,blue-green-yellow-red
%  [34,35,36]:AFM hot,black-red-yellow-white
%  [ 3, 3, 3]:gray
%
% version:201800402

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colormap formulae:                                       %
%                                                          %
% Copyright (C) 2018 Masayuki Tanaka. All rights reserved. %
%                    mtanaka@sc.e.titech.ac.jp             %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ids = type2ids( type )

 if( ~exist('type','var') || isempty(type) )
  fprintf('[ 7, 5,15]:traditional,pm3d,black-blue-red-yellow\n');
  fprintf('[ 3,11, 6]:green-red-violet\n');
  fprintf('[23,28, 3]:ocean,green-blue-white\n');
  fprintf('[21,22,23]:hot,black-red-yellow-white\n');
  fprintf('[30,31,32]:color printable on gray,black-blue-violet-yellow-white\n');
  fprintf('[33,13,10]:rainbow,blue-green-yellow-red\n');
  fprintf('[34,35,36]:AFM hot,black-red-yellow-white\n');
  fprintf('[ 3, 3, 3]:gray\n');
  

 else
 
  if( strcmpi( type, 'traditional' ) || strcmpi( type, 'pm3d' ) || strcmpi( type, 'black-blue-red-yellow' ) )
   ids = [7,5,15];
  elseif( strcmpi( type, 'green-red-violet' ) )
   ids = [3,11,6];
  elseif( strcmpi( type, 'ocean' ) || strcmpi( type, 'green-blue-white' ) )
   ids = [23,28,3];
  elseif( strcmpi( type, 'hot' ) || strcmpi( type, 'black-red-yellow-white' ) )
   ids = [21,22,23];
  elseif( strcmpi(type, 'color printable on gray') || strcmpi( type, 'black-blue-violet-yellow-white' ) )
   ids = [30,31,32];
  elseif( strcmpi(type, 'rainbow') || strcmpi( type, 'blue-green-yellow-red' ) )
   ids = [33,13,10];
  elseif( strcmpi(type, 'AFM hot') || strcmpi( type, 'black-red-yellow-white' ) )
   ids = [34,35,36];
  elseif( strcmpi(type, 'gray') )
   ids = [3,3,3];
  else
   error( 'Unknown type' );
  end
 
 end
end
