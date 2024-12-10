function [output] = vec2list(A,times)
if numel(A)~=numel(times)
    error('Size of Lists do not match.')
end
if ~isrow(A)
    A = A';
end
if ~isrow(times)
    times = times';
end
B = find(A);
x = diff(B)==1;
f = find([false,x]~=[x,false]);
g = find(f(2:2:end)-f(1:2:end-1)>=1);
startI = B(f(2*g-1));
endI = zeros(size(startI));
for i = 1:numel(startI),
    endP = find(A(startI(i):end)==0,1,'first');
    if ~isempty(endP)
        endI(i) = startI(i)+endP-2;
    else
        endI(i) = length(A);
    end
end
output = [times(startI)',times(endI)'];

end