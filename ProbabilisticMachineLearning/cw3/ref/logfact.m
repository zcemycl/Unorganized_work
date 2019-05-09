function lf = logfact(num)
dim = size(num);
if dim(1)>dim(2)
    dim=dim;
elseif dim(2)>dim(1)
    dim=dim';
end
lf = [];
for i = 1:dim(1)
lf(i) = sum(log(linspace(1,num(i),num(i))));
end
end