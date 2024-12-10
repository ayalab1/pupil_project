function [xvect,output] = wb_list2vec(A,times)
output = zeros(size(times));
xvect = [];

for i=1:length(A(:,1)),
    sp = find(times>=A(i,1),1,'first');
    ep = find(times<=A(i,2),1,'last');
    if ~isempty(sp)
        output(sp:ep) = 1;
        xvect(i,1) = times(sp);
        xvect(i,2) = times(ep);
    end
end