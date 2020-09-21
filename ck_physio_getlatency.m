
function Latency = ck_physio_getlatency(sigZ,tax,thr,INTERVAL,DUR)

% function Latency = ck_physio_getlatency(sigZ,tax,thr,INTERVAL,DUR)
%
% function to estimate response threshold as well as latency using the SD method
% INTERVAL: Duration of window for which resp > thr
% DUR: maximal duration of stimulus interval

if nargin < 5
  DUR = 800;
end
if nargin < 4
  INTERVAL = 20; % msec
end
if nargin < 3
  thr = [2:0.5:3];
end

rate = 1000/median(diff(tax));

for T=1:length(thr)
  X = (sigZ>thr(T));
  X2 = ck_helper_smoothsig(X,rate,INTERVAL,2);
  % assume time interval for search
  tmax = max(tax,DUR);
  J = find( (tax>=-100).*(tax<=tmax));
  X2 = (X2(J)>0.92)+(X2(J)<-0.92);  
  t = find(X2);
  if ~isempty(t)
    t = min(t);
    Latency(T,:) = [t,1];
  else
    Latency(T,:) = [0,0];
  end
end

Latency(:,1) = Latency(:,1) * (1000/rate) - 100;
Latency(:,1) = Latency(:,1) - INTERVAL/2;
return;

