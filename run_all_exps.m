
List = {'I05Da','I05Dc','I05Dd','I05Df','I05Dh','I05Di','I05Dj','I05Dl','I05Dm','I05Em','I05En','I05Ep','I05Eq','I05Er','I05Es','I05Ew','I05Ex'};

for k=1:length(List)
  ck_physio_load(List{k},[]);
  ck_physio_anasig(List{k});
end
