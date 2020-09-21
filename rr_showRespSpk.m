
function rr_showRespSpk(Ses,grp,N_Sites)

% function to load & produce PSTH figure of recording sites with depth

eval(Resp.Ses);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation
Colorvec = 'rgbcmykr--g--b--c--m--y--k--rgbcmykrgbcmykrgbcmykrgbcmyk';
N_el = length(Resp.Spk);
N_stim = length(Resp.Lfp{1}.Evp);
triallength=
binsize=
fs=      %sampling rate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computing

lastBin=binsize*ceil((triallength-1)*(1000/(fs*binsize)));
limits=0:binsize:lastBin;
x=(mod(times-1,triallength)+1)*(1000/fs);
r=(histc(x,limits)*1000)/(ntrials*binsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting

axes(h);
ph=bar(limits(1:end-1),r(1:end-1),'hist');
set(ph,'edgecolor',Colorvec,'facecolor',Colorvec);