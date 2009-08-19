function vcombs=combs(varargin)
%function vcombs=combs(varargin)
%
%returns matrix of all combinations of input values
%size of output is prod(l1...ln) x n
%where vargin contains n arrays with lengths l1, l2, ...
%if there are 2 arguments and the second is scalar p, output is same as for p repetitions of 1st argument

dims=numel(varargin);
if dims==2
    for n=2:varargin{2}
        varargin(n)=varargin(1);
    end
end
dims=numel(varargin);

for d=1:dims
	nvals(d)=numel(varargin{d});
end
ncombs=prod(nvals);
denom=[ flipdim(cumprod(flipdim(nvals(2:end),2)),2) 1];

ind=1:ncombs;
vcombs=zeros(ncombs,dims);
for d=1:dims
   vcombs(:,d)= varargin{d}(floor(1+mod((ind-1)/denom(d),nvals(d))));
end
