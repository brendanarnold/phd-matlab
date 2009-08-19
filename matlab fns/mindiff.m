function r=mindiff(a,b)
%returns matrix, the size of a, which contains the difference |a_i-b_j|
%where b_j is element of b that is closest to a_i

m1=repmat(a(1:end)',1,numel(b));
m2=repmat(b(1:end),numel(a),1);
r=min(abs(m1-m2),[],2);
r=reshape(r,size(a));