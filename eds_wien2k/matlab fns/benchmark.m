t1=clock;
a=rand(1000,1000);
b=inv(a);
t2=clock;
disp(etime(t2,t1));
