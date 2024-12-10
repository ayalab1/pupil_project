function [lo,hi]= findcontiguous(x, order)
%function [lo,hi]= findcontiguous(x, order)
%  
%   Find the beginnings (lo) and ends (hi) of contiguous values in x.
%   E.g. x= [1 2 5 6 7 8 10] results in
%       lo= [1 5 10]
%       hi= [2 8 10]

lo= []; hi=[];
if ~isempty(x); 
    if nargin<2; order= 1; end

    if ~any(size(x)==1); error('function only works for vectors'); end
    row= 0;
    if size(x,1)==1; 
        row= 1;
        x= x';
    else
        row= 0;
    end

    dind= find(diff(x)>order);
    nind= length(x);
    hi= [x(dind); x(nind)]+order-1;
    lo= x([1; dind+1]);

    if(row)
        lo= lo';
        hi= hi';
    end
end
