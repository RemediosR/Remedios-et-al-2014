
X = zeros(1,20833);
Y = zeros(1,80);
Y(41:50) = +1;
Y(30:40) = -1;
Y = (conv(Y,[1,2,3,4,5,5,4,3,2,1]));

dt = ceil(20833/1000);
c = 1;
while c<(length(X)-100)
  del = rand*200;
  del = round(del*dt);
  X(c+del:c+del+(length(Y)-1)) = Y;
  c=c+del;
end

save('D:\MpiPhysiology\testSpk.mat','X');
