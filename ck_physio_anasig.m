function ck_physio_anasig(Ses,grp,site)

% function ck_physio_anasig(Ses,stimulus,site)
%
% function to analyze the band limited signals for an experiment or a group
% analysis LFP, AMUA, SPK and Eye movements if present
% 

% Stimulus set dependent parameters 

ARG.PAR.No_time_adjust = {'flash','whisker','somato'};
ARG.PAR.Do_revcor = {'fastmonk','fastnat','fastfaces','voices_faceslong','fastsounds2','faceview','faceview2'};
ARG.PAR.Specialized = {'monkadapt'};

% This is a list of paradigms for which pure visual correction is used.
ARG.PAR.Timing.visual = {'fastfaces','faces','faceview','faceview2'};
% a list of paradigms where auditory / or visual correction is used
% dependent on stimulus
ARG.PAR.Timing.mixed = {'tpoaudvis','audvis','audvisIE','audvismatch','audvisnoise','audvis2','audviscp','voices_faces','voices_faceslong'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for signals
ARG.Lfp.unit = 'sd'; % units for spectrum
ARG.Lfp.nboot = 1000;
ARG.Lfp.spec.nfft = 256; % nfft for spectrum
ARG.Lfp.spec.olp = 195; % overlap 
ARG.Lfp.spec.winl = 200;  % window length in ms.
ARG.wavelet.Freqs = [5,10,20,40,60,80,120];
ARG.wavelet.nco = 4.5; % number of side lobes
ARG.Lfp.pval = 0.05; % p-value for bootstrapping

ARG.Amua.unit = 'sd'; % units for Amua data
ARG.Amua.nboot = 1000;
ARG.Amua.pval = 0.05;

ARG.Spk.unit = 'imp';
ARG.Spk.bin = 5; % temporal binning of spikes in ms
ARG.Spk.psth_smooth = 20; % msec kernel for smoothing
ARG.Spk.hr_latency = 1; % compute latency at 2ms resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other parameters
ARG.etc.time_adjust = 1; % timing adjustment for the different stimuli (options 0/1/2)
ARG.etc.get_eye = 1; % get eye movement data

ARG.etc.bad_trial_thr_lfp = 4; % threshold for kicking trials based on LFP variance
ARG.etc.bad_trial_thr_amua = 4; % threshold for kicking trials based on AMUA variance
% LARGE VALUES MEAN LITTLE THRESHOLDING. If zero, no thresholding is done

ARG.Stim_int = [50,200]; % stimulus interval to compute response amplitude and responsiveness

% some hack for specgram function
global spec_flag
if exist('spectrogram')==2
  spec_flag = 1;
else spec_flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% files to analyze
INFO = ck_physioprepare(Ses);
eval(Ses);
if nargin ==1
  % do all groups and sites
  Do_groups = [1:INFO.NF];
  for g=1:length(Do_groups)
    Do_sites{g} = [1:INFO.NSITES(g)];
  end
elseif nargin==2
  % do all sites for this group
  tmp = strmatch(grp,INFO.fields,'exact');
  if isempty(tmp)
    fprintf('%s does not seem to be a valid groupd\n',grp);
    return;
  end    
  Do_groups = tmp;
  Do_sites{1} = [1:INFO.NSITES(tmp)];
elseif nargin==3
    tmp = strmatch(grp,INFO.fields,'exact');
  if isempty(tmp)
    fprintf('%s does not seem to be a valid groupd\n',grp);
    return;
  end    
  Do_groups = tmp;
  Do_sites{1} = site;
end
PATHNAME = SYSP.DataNeuro;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  loop groups and sites
% This is kind of a mess since several parameters are reset for specific
% paradigms, different evaluation functions called, paradigms omitted
% etc...
fprintf('analyzing %s...\n',Ses);
for G=1:length(Do_groups)
  grpname = INFO.fields{Do_groups(G)};
  if strfind(grpname,'longnat')
    fprintf('Skipping %s\n',grpname);
    continue;
  end
  % make sure anasig is never run on longnat
  for S=1:length(Do_sites{G})
    ARGL = ARG;
    % some special cases
    if sum(strcmp(grpname,{'fastMonk','fastNat','fastFaces','FastSounds2'}))
      ARGL.etc.bad_trial_thr_lfp = 7; 
      ARGL.etc.bad_trial_thr_amua = 7; 
    end
    % stop time adjustment for somatosensory paradigms
    if sum(strcmp(grpname,ARG.PAR.No_time_adjust))
      ARGL.etc.time_adjust = 0;
      if sum(strcmp(grpname,{'flash'}))
        ARGL.etc.time_adjust = 2; % special time adjustment for always AV paradigms
      end
    end    
    fprintf('group   %s \t site %d \t',grpname,Do_sites{G}(S));
    Filename =  sprintf('%s/%s_%s_site%02d.mat',PATHNAME,SYSP.dirname,grpname,Do_sites{G}(S));
    DD = load(Filename);
    % adjust paradigm for revcor experiments
    if isfield(DD.Data{1}.info.BEH,'RevCor')
      if sum(strcmp(lower(grpname),ARG.PAR.Do_revcor))
        DD.Data = ck_physio_fastParad(DD.Data,grpname);
        ARGL.etc.time_adjust = 1;
      end
    end
    if sum(strcmp(lower(grpname),ARG.PAR.Timing.visual))
      ARGL.timing = 2;
    elseif sum(strcmp(lower(grpname),ARG.PAR.Timing.mixed))
      ARGL.timing = 3;
    else
      ARGL.timing = 1;
    end

    % get the responses - for some paradigms we have to defer to specialized
    % subfunctions
    
    if sum(strcmp(lower(grpname),ARG.PAR.Specialized))
      % call the respective specialized functions
      if strfind(lower(grpname),'adapt')
        funname = 'SubGetAdapt';
      end
      eval(['Resp = ' funname '(DD.Data,ARGL);']);
    else
      % default
      Resp = subGetResponse(DD.Data,ARGL);
    end
    Resp.grp = grpname;
    Resp.site = Do_sites{G}(S);
    Resp.ARG = ARGL;
    Resp.Ses = DD.Data{1}.Ses;
    Resp.info = DD.Data{1}.info;
    % save result
    sname = sprintf('%s/%s_%s_ANALYZED_site%02d.mat',PATHNAME,SYSP.dirname,grpname,Do_sites{G}(S));
    save(sname,'Resp');
    fprintf('\n');
  end
end


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Resp = subGetResponse(Data,ARG)

% comptue timingn adjustments for sound presentation
% if not set, trivial offset is returned
[time_adjust,computed_adjust] = ck_physio_gettrialtiming(Data,ARG);
adjust.time_adjust = time_adjust;
adjust.computed_adjust = computed_adjust;
Resp.adjust = adjust;

if ARG.etc.get_eye 
  Resp.Eye = subGetEye(Data,ARG,time_adjust);
end

if ARG.etc.bad_trial_thr_lfp
  % compute variability of signal to remove (bad) trials
  Resp.data_quality = ck_physio_DataQ(Data,ARG,time_adjust);
else
  Resp.data_quality.mask =[];
  Resp.data_quality.num_bad = zeros(length(Data),length(Data{1}.info.BEH.conds));
end

% compute the Lfp, Amua and SPK responses
if ~isfield(Data{1},'Lfp')
  fprintf('Lfp not found');
else
  Resp.Lfp = subGetLfp(Data,ARG,adjust,Resp.data_quality.mask);
end

if ~isfield(Data{1},'Amua')
  fprintf('Amua not found'); 
else
  Resp.Amua = subGetAmua(Data,ARG,adjust,Resp.data_quality.mask);
end
if ~isfield(Data{1},'Spk')
  fprintf('Spks not found'); 
  for k=1:length(Resp.Lfp)
    Resp.Spk{k}.nunits = 0;
  end
else
  Resp.Spk = subGetSpk(Data,ARG,adjust,Resp.data_quality.mask);
  if ARG.Spk.hr_latency
    tmplat = subGetSpkLat(Data,ARG,adjust,Resp.data_quality.mask);
    for ee=1:length(tmplat)
      if isfield(tmplat{ee},'sigrespZ')
        Resp.Spk{ee}.Spk_lat = tmplat{ee}.sigrespZ;
      end
    end
  end
end

fprintf('   Disc. trials/elect: ');
for E=1:length(Data)
  fprintf('%d  ',sum(Resp.data_quality.num_bad(E,:)));
end

return;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Lfp = subGetLfp(Data,ARG,adjust,MASK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_adjust = adjust.time_adjust;
computed_adjust = adjust.computed_adjust;

if isempty(MASK)
  do_mask = 0;
else
  do_mask = 1;
end

global spec_flag
fprintf('Lfp\t');
% get lfp response
Nel = length(Data);
rate = Data{1}.info.ARG.Lfp.rate;
for E=1:Nel
  dat = Data{E}.Lfp;
  model = Data{E}.Lfp_model;
  Nstim = length(dat);
  t_ax = Data{E}.Lfp_taxis;
  T_B = find(Data{Nel}.Lfp_taxis<=-30);
  T_S = find(model==1); % all stimulus duration

  % compute the evoked potential
  for S=1:Nstim 
    if isempty(dat{S})
      continue;
    end
    tmp = dat{S};
    Lfp{E}.Obs_id{S} = single(Data{E}.Lfp_obsid{S});
    % time-adjust individual trials
    offset = round(time_adjust{S}/1000*rate);
    tmp = subDataOffsetCorrect(tmp,offset);
    % write back for later analysis
    dat{S} = tmp;
    if do_mask
      % KICK TRIALS
      good = find(MASK{E,S}==0);
      tmp = tmp(good,:);
    end
    if min(size(tmp))==0
      continue;
    end
    tmp = subDataPrep(tmp,ARG,T_B,'Evp');
    evp = mean(tmp,1);
    Lfp{E}.Evp{S} = evp;
    if ARG.Lfp.nboot
      t1=bootstrp(ARG.Lfp.nboot,'mean',tmp);
      Ci_lo = prctile(t1,ARG.Lfp.pval*100);
      Ci_up = prctile(t1,100-ARG.Lfp.pval*100);
      Lfp{E}.Evp_ci{S}(1,:) = Ci_lo;
      Lfp{E}.Evp_ci{S}(2,:) = Ci_up;
    end
    
    % some response measure for tuning curves
    T_J = find( ( t_ax>=ARG.Stim_int(1)).*( t_ax<=ARG.Stim_int(2)));
    Lfp{E}.tuning_evp(S) = mean(mean(abs(tmp(:,T_J))),2);
    Lfp{E}.tuning_evps(S) = std(mean(abs(tmp(:,T_J)),2))/sqrt(size(tmp,1));
    
    % test for responsiveness
    T_B2 = find( (t_ax <= -50).*(t_ax>=-300));
    base = mean(abs(tmp(:,T_B2)),2);
    stim = mean(abs(tmp(:,T_J)),2);
    warning off
    [h,p] = ttest(base,stim);
    warning on
    Lfp{E}.Evp_sigresp(S) = p/2;
    % z-score signal for latency analysis
    evp2 = subDataPrep(mean(tmp,1),ARG,T_B,'Evp');
    Lfp{E}.Evp_sigrespZ(S,:,:) = ck_physio_getlatency(evp2,t_ax);
    
  end
  Lfp{E}.Evp_model = model;
  Lfp{E}.Evp_taxis = Data{E}.Lfp_taxis;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute spectrogram for each trial
  % one dummy run to get matrix size
  dat = Data{E}.Lfp;
  % adjust parameters to the current sampling rate
  swinl = round(ARG.Lfp.spec.winl/(1000/rate));
  solp = round(ARG.Lfp.spec.olp/(1000/rate));
  
  for S=1:Nstim
    if isempty(dat{S})
      continue;
    end
    tmp = dat{S}(1,:);
    if spec_flag
      [spec,fax,tax] = spectrogram(tmp,swinl,solp,ARG.Lfp.spec.nfft,rate);
    else
      [spec,fax,tax] = specgram(tmp,ARG.Lfp.spec.nfft,rate,swinl,solp);
    end
    tax = tax*1000 + Lfp{E}.Evp_taxis(1);
    tax = tax + (swinl);
    Lfp{E}.spec_t{S} = tax;
    T_Bs = find(tax<+10);
    if isempty(T_Bs)
      T_Bs = find(tax<+30);
    end
    Specs = zeros(size(dat{S},1),size(spec,1),size(spec,2));
    for trial = 1:size(dat{S},1)
      tmp = dat{S}(trial,:);
      if spec_flag
        [spec] = spectrogram(tmp,swinl,solp,ARG.Lfp.spec.nfft,rate);
      else
        [spec] = specgram(tmp,ARG.Lfp.spec.nfft,rate,swinl,solp);
      end
      spec = abs(spec);
      % change units
      spec =  subDataPrep(spec,ARG,T_Bs,'Lfp');
      Specs(trial,:,:) = spec;
    end
    if do_mask
      % KICK TRIALS
      good = find(MASK{E,S}==0);
      Specs = Specs(good,:,:);
    end
    if min(size(Specs))==0
      continue;
    end

    
    model2 = zeros(size(tax));
    T_stim = find( (tax>(T_S(1)/rate)).*(tax<T_S(end)/rate));
    model2(T_stim) =1 ;
    % keep only frequencies up to 120 Hz
    tmp = find(fax<=120);
    fax = fax(tmp);
    Specs = Specs(:,tmp,:);
    Lfp{E}.Spec{S} = squeeze(mean(Specs,1));
    Lfp{E}.spec_t{S} = tax;
    Lfp{E}.spec_f{S} = fax;
    Lfp{E}.spec_model = model2;
  end % stim
  
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute spectrogram for each trial using wavelets
  % generate wavelets
  Freq = ARG.wavelet.Freqs;
  for f=1:length(Freq)
    s = ARG.wavelet.nco/(9*Freq(f));
    X = [-20:1/rate:20];
    GG = exp( -(X.^2/(2*s^2))).*exp(i*2*pi*Freq(f)*X);
    % truncate
    zx{f} = GG-mean(GG);
    GG = GG(find(abs(GG)>0.01));
    % zero shift, not so important
    GG = GG-mean(real(GG));
    Wave{f} = GG;
  end
  
  % compute wavelet filtered data
  t_ax = Data{E}.Lfp_taxis;
  model = Data{E}.Lfp_model;
  T_B = find((t_ax<=-30));
  T_S = find(model);
  Lfp{E}.wave_tax = Data{E}.Lfp_taxis;
  
  for f=1:length(Wave)
    n1 = floor(length(Wave{f})/2);
    n2 = ceil(length(Wave{f})/2);
    for S=1:Nstim
      if isempty(dat{S})
        continue;
      end
      X = dat{S};
      ntrial = size(X,1);
      Dfilt = zeros(ntrial,size(X,2));
      for t=1:ntrial
        Y2 = X(t,:);
        % extend window to prevent edge-artifacts at low freqs
        % might be better than zero padding
        Y2 = [Y2(end:-1:1) Y2 Y2(end:-1:1)];
        tmp = conv(Y2,Wave{f});
        tmp = tmp(n1:end-n2);
        lle = size(X,2);
        tmp = tmp([lle+1:2*lle]);
        Dfilt(t,:) = tmp;
      end
      if do_mask
        % KICK TRIALS
        good = find(MASK{E,S}==0);
        Dfilt = Dfilt(good,:);
      end
      if min(size(Dfilt))==0
        continue;
      end

      
      Dfiltamp = subDataPrep(abs(Dfilt),ARG,T_B,'wave');
      Lfp{E}.wave{S,f} = mean(Dfiltamp,1);
      % phase consistency
      % Lfp{E}.wave_phasecons{S,f} = mean(exp(i*angle(Dfilt)),1);
      % wave amplitude
      J = find( (t_ax>ARG.Stim_int(1)).*(t_ax<=ARG.Stim_int(2)));
      act = mean(Dfiltamp(:,J),2);
      Lfp{E}.wave_tunign_m(S,f) = mean(act);
      Lfp{E}.wave_tunign_s(S,f) = std(act)/sqrt(size(act,1));
    end % S
  end % f
  
  Lfp{E}.elec = Data{E}.elec;
  Lfp{E}.grpname = Data{E}.info.grpname;
  Lfp{E}.site = Data{E}.info.site;  
end


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Amua = subGetAmua(Data,ARG,adjust,MASK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
time_adjust = adjust.time_adjust;
computed_adjust = adjust.computed_adjust;
if isempty(MASK)
  do_mask = 0;
else
  do_mask = 1;
end

% get amua response
fprintf('Amua\t');
Nel = length(Data);
rate = Data{1}.info.ARG.Amua.rate;
for E=1:Nel
  dat = Data{E}.Amua;
  model = Data{E}.Amua_model;
  Nstim = length(dat);
  t_ax = Data{E}.Amua_taxis;
  T_B = find(t_ax<=-30);
  T_S = find(model==1);
  for S=1:Nstim 
    if isempty(dat{S})
      continue;
    end
    % time-adjust individual trials
    tmp = dat{S};
    offset = round(time_adjust{S}/1000*rate);
    tmp = subDataOffsetCorrect(tmp,offset);    
    dat{S} = tmp;
    if do_mask
      good = find(MASK{E,S}==0);
      tmp = tmp(good,:);
    end
    if min(size(tmp))==0
      continue;
    end

    tmp = subDataPrep(tmp,ARG,T_B,'Amua');  
    Amua{E}.Evp{S} = mean(tmp,1);
    if ARG.Amua.nboot
      t1=bootstrp(ARG.Amua.nboot,'mean',tmp);
      Ci_lo = prctile(t1,ARG.Amua.pval*100);
      Ci_up = prctile(t1,100-ARG.Amua.pval*100);
      Amua{E}.Evp_ci{S}(1,:) = Ci_lo;
      Amua{E}.Evp_ci{S}(2,:) = Ci_up;
    end
    Amua{E}.t_axis = Data{E}.Amua_taxis;
    
    % some stuff for tuning
    T_J = find( (Data{E}.Amua_taxis>=ARG.Stim_int(1)).*(Data{E}.Amua_taxis<=ARG.Stim_int(2)));
    Amua{E}.tuning_m(S) = mean(mean(abs(tmp(:,T_J)),2));
    Amua{E}.tuning_s(S) = std(mean(abs(tmp(:,T_J)),2),[],1)/sqrt(size(tmp,1));

    T_B2 = find( (t_ax <= -50).*(t_ax>=-300));
    base = mean( abs(tmp(:,T_B2)) ,2);
    stim = mean( abs(tmp(:,T_J)) , 2);
    warning off
    [h,p] = ttest(base,stim);
    warning on
    Amua{E}.Evp_sigresp(S) = p/2;    
    evp2 = subDataPrep(mean(tmp,1),ARG,T_B,'Amua');
    Amua{E}.Evp_sigrespZ(S,:,:) = ck_physio_getlatency(evp2,t_ax);
  end
 
    
  Amua{E}.elec = Data{E}.elec;
  Amua{E}.grpname = Data{E}.info.grpname;
  Amua{E}.site = Data{E}.info.site;
  Amua{E}.model = model;
end
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SPK = subGetSpk(Data,ARG,adjust,MASK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
time_adjust = adjust.time_adjust;
computed_adjust = adjust.computed_adjust;

% Spk field: Timestamp, PCA1, PCA2, PCA3, nonlienar energy
% Spk{stim}{unit}{trial}
% Spk_detail: Channel, Unit, Timestamp, PCA1, PCA2, PCA3, nonlienar energy, waveform, spk_file number
if isempty(MASK)
  do_mask = 0;
else
  do_mask = 1;
end

fprintf('Spk\t');
Nel = length(Data);
bin = ARG.Spk.bin;
for E=1:Nel
  dat = Data{E}.Spk;
  if ~isempty(dat)
    Nstim = length(dat);
    trial_dur = ceil( abs(diff(Data{E}.Lfp_taxis([1,end]))) / bin)+5;
    Nunits = length(dat{1});
    for U=1:Nunits
      for S=1:Nstim
        if isempty(dat{S})
          continue;
        end
        tmp = dat{S}{U};
        Ntrial = length(tmp);
        Matrix = zeros(Ntrial,trial_dur);
        for t=1:Ntrial
          if ~isempty(tmp{t})
            spk_time = tmp{t}(:,1)*1000; % units of ms
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % time-adjust spike timing !!!!!!!!!!!!!!!!
            % this  might create spikes before trial onset, remove these
            spk_time = spk_time - time_adjust{S}(t);  
            spk_time = ceil(spk_time / bin); % binning
            spk_time = spk_time(find(spk_time>0));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for l=1:length(spk_time)
              Matrix(t,spk_time(l)) = Matrix(t,spk_time(l)) + 1;
            end
          end % if
        end % trial
        Matrix = Matrix * (1000/bin); % THIS ONLY FOR UNITS OF SPK/sec
        if do_mask
          good = find(MASK{E,S}==0);
          Matrix = Matrix(good,:);
        end
        if min(size(Matrix))==0
          continue;
        end

        SPK{E}.resp{U}{S} = Matrix;
        SPK{E}.Obs_id{U}{S} = single(Data{E}.Lfp_obsid{S}(good,:));
        width = ARG.Spk.psth_smooth/ARG.Spk.bin;
        width = ceil(width*1.8);
        Matrix = sub_smooth_spk(Matrix,width);
        
        % confidence interval for mean       
        if Ntrial ==1
          Matrix(2,:) = Matrix(1,:);
        end
        if ARG.Amua.nboot
          t1=bootstrp(ARG.Amua.nboot,@mean,Matrix);
          Ci_lo = prctile(t1,ARG.Amua.pval*100);
          Ci_up = prctile(t1,100-ARG.Amua.pval*100);
          SPK{E}.mean_ci(U,1,:) = Ci_lo;
          SPK{E}.mean_ci(U,2,:) = Ci_up;
        end
        SPK{E}.mean_resp{U}(S,:) = mean(Matrix,1);
        
        model = zeros(1,size(SPK{E}.mean_resp{U},2));
        t_axis = [1:length(model)] -1 + round(Data{E}.Spk_taxisoffset/bin);
        t_axis = t_axis*bin;
        stim_dur = (Data{E}.info.Cut.dur_stim);
        model(find( (t_axis>0).*(t_axis <= stim_dur))) = 1;        
        SPK{E}.model = model;
        SPK{E}.t_axis = t_axis;
        
        % stuff for tuning
        T_J = find( (t_axis>=ARG.Stim_int(1)).*(t_axis<=ARG.Stim_int(2)));
        resp = mean(Matrix(:,T_J),2);
        SPK{E}.resp_m(U,S) = mean(resp);
        SPK{E}.resp_s(U,S) = std(resp)/sqrt(size(resp,1));
           
        % t-test for response
        % use baseline vs. first 400ms of stimulus
        T_B2 = find( (t_axis <= -50).*(t_axis>=-300));
        val1 = sum(Matrix(:,T_B2),2);
        val2 = sum(Matrix(:,T_J),2);
        % now we have to make sure that it works for all paradigms
        if length(val1)>1
          if var(val1)==0
            val1(1) = val1(2)+0.01;
          end
          if var(val2)==0
            val2(1) = val2(2)+0.01;
          end
          warning off
          [h,p] = ttest2(val1,val2);
          warning on
          SPK{E}.sig_resp(U,S) = p/2;
        else
          SPK{E}.sig_resp(U,S) = 1;
        end
        T_B =  find(t_axis<=0);
        evp2 = subDataPrep(mean(Matrix,1),ARG,T_B,'SPK');        
        SPK{E}.sigrespZ(U,S,:,:) = ck_physio_getlatency(evp2,t_axis);
      end % S
      % store spike waveform
      Wave = Data{E}.Spk_detail{U}(:,8:end-1);
      if ~isempty(Wave)
        xmean = mean(Wave,1);
        sem = std(Wave)/sqrt(size(Wave,1));
        crit = tinv((1 - 0.05 / 2), (size(Wave,1)-1)) .* sem;
        SPK{E}.wave{U}(1,:) = xmean;
        SPK{E}.wave{U}(2,:) = xmean - crit;
        SPK{E}.wave{U}(3,:) = xmean + crit;
        SPK{E}.wave{U}(4,:) = std(Wave);
        % signal to noise ratio of wave form
        b = std(xmean([1:5,end-1:end]));
        p = max(abs(xmean));
        SPK{E}.wave_snr(U) = p/b;
        % ISI histogram with 0.5 ms bin
        isi = Data{E}.Spk_detail{U}(:,3);
        isi = diff(isi*1000);
        isi = hist(isi,[0:0.5:400]);
        SPK{E}.wave_ISI{U} = isi/sum(isi);
      else
        SPK{E}.wave{U} = [];
        SPK{E}.wave_snr(U) = 0;
        SPK{E}.wave_ISI{U} =[];
      end
   end % units
   SPK{E}.nunits =  Nunits;
  else
    SPK{E}.nunits = 0;
  end % if
end



return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LAT = subGetSpkLat(Data,ARG,adjust,MASK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute latency for 1ms resolution
time_adjust = adjust.time_adjust;
computed_adjust = adjust.computed_adjust;
if isempty(MASK)
  do_mask = 0;
else
  do_mask = 1;
end

Nel = length(Data);
bin = 2; % 1ms resolution for good spike latency
psth_smooth = 10; %ms 

for E=1:Nel
  dat = Data{E}.Spk;
  if ~isempty(dat)
    Nstim = length(dat);
    trial_dur = ceil( abs(diff(Data{E}.Lfp_taxis([1,end]))) / bin)+5;
    Nunits = length(dat{1});
    for U=1:Nunits
      for S=1:Nstim
        if isempty(dat{S})
          continue;
        end
        tmp = dat{S}{U};
        Ntrial = length(tmp);
        Matrix = zeros(Ntrial,trial_dur);
        for t=1:Ntrial
          if ~isempty(tmp{t})
            spk_time = tmp{t}(:,1)*1000; % units of ms
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % time-adjust spike timing !!!!!!!!!!!!!!!!
            % this  might create spikes before trial onset, remove these
            spk_time = spk_time - time_adjust{S}(t);  
            spk_time = ceil(spk_time / bin); % binning
            spk_time = spk_time(find(spk_time>0));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for l=1:length(spk_time)
              Matrix(t,spk_time(l)) = Matrix(t,spk_time(l)) + 1;
            end
          end % if
        end % trial
        Matrix = Matrix * (1000/bin); % THIS ONLY FOR UNITS OF SPK/sec
        if do_mask
          good = find(MASK{E,S}==0);
          Matrix = Matrix(good,:);
        end
        if min(size(Matrix))==0
          continue;
        end

        SPK{E}.resp{U}{S} = Matrix;
        SPK{E}.Obs_id{U}{S} = single(Data{E}.Lfp_obsid{S}(good,:));
        width = psth_smooth/ARG.Spk.bin;
        width = ceil(width*1.8);
        Matrix = sub_smooth_spk(Matrix,width);
        
        % confidence interval for mean       
        if Ntrial ==1
          Matrix(2,:) = Matrix(1,:);
        end
        
        model = zeros(1,size(Matrix,2));
        t_axis = [1:length(model)] -1 + round(Data{E}.Spk_taxisoffset/bin);
        t_axis = t_axis*bin;
        T_B =  find(t_axis<=0);
        evp2 = subDataPrep(mean(Matrix,1),ARG,T_B,'SPK');        
        LAT{E}.sigrespZ(U,S,:,:) = ck_physio_getlatencyHR(evp2,t_axis);
      end % S
    end % units
    LAT{E}.nunits =  Nunits;
  else
    LAT{E}.nunits = 0;
  end % if
end



return;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the eye movement data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Eye = subGetEye(Data,ARG,time_adjust);

Nstim = length(Data{1}.Lfp);
EM = Data{1}.info.BEH.em;
tax = Data{1}.Lfp_taxis;
Ton = Data{1}.info.Cut.Time_on;
dur = ceil((tax(end)-tax(1))/5);

if isempty(EM{1})
  fprintf('EYE MOVEMNTS NOT FOUND\n');
  Eye =[];
  
  return;
end

for S=1:Nstim
  oid = Data{1}.Lfp_obsid{S};
  trials = oid(:,1) + 1; % +1 for true_obs
  for k=1:length(trials)
    X = EM{trials(k)}{2};
    Y = EM{trials(k)}{3};
    t_eye = ([1:length(X)]-1) * EM{trials(k)}{1};
    onset = Ton(trials(k));
    [val,ind] = min(abs(t_eye-(onset+tax(1))));
    t_eye = t_eye(ind:end);
    X = X(ind:end); Y = Y(ind:end);
    t_eye = t_eye - t_eye(1)+tax(1);
  
    if length(t_eye)>=dur
      t_eye = t_eye(1:dur);
      X = X(1:dur); Y = Y(1:dur);
    else
      % this should happen for only a very few trials where the state system dropped some datapoints
      dt = (t_eye(2)-t_eye(1));
      t_eye = [t_eye(1):dt:t_eye(1)+(dur-1)*dt];
      X = zeros(1,dur); Y = zeros(1,dur);      
    end
    Eye{S}.dat(k,1,:) = single(X);
    Eye{S}.dat(k,2,:) = single(Y);
    Eye{S}.tax = t_eye;
    Eye{S}.calib = Data{1}.info.BEH.eyecalib;
     
  end
end

  


return;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dat = subDataPrep(dat,ARG,blank,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch mode
  case {'Evp'} %%%%%%%%%%%%%%%%%
    % remove mean in baseline
    dat = dat - repmat(mean(dat(:,blank),2),[1,size(dat,2)]);
    % z-score units
    dat = dat./repmat(std(dat(:,blank),[],2),[1,size(dat,2)]);
  case {'wave'} %%%%%%%%%%%%%%%%%
    base = repmat(mean(dat(:,blank),2),[1,size(dat,2)]);
    switch ARG.Lfp.unit
      case {'sd'}
       dat = dat - base;
       dat = dat./repmat(std(dat(:,blank),[],2),[1,size(dat,2)]);
     case {'pc'}
       dat = (dat - base)./base*100;
    end
  case {'Lfp'}  %%%%%%%%%%%%%%%%%
   % normalize the spectrogram
   switch ARG.Lfp.unit
     case {'sd'}
       dat = dat - repmat(mean(dat(:,blank),2),[1,size(dat,2)]);
       dat = dat./repmat(std(dat(:,blank),[],2),[1,size(dat,2)]);
     case {'pc'}
       tmp  =repmat(mean(dat(:,blank),2),[1,size(dat,2)]);
       dat = (dat-tmp)./tmp*100;       
   end
  case {'Amua'} %%%%%%%%%%%%%%%%%
    switch ARG.Amua.unit
      case {'sd'}
        % remove mean in baseline
        dat = dat - repmat(mean(dat(:,blank),2),[1,size(dat,2)]);
        % z-score units
        dat = dat./repmat(std(dat(:,blank),[],2),[1,size(dat,2)]);
      case {'pc'}
        base = repmat(mean(dat(:,blank),2),[1,size(dat,2)]);
        dat = (dat-base)./base*100;
    end
    
  case {'SPK'}
    % SD units    
    dat = dat - repmat(mean(dat(:,blank),2),[1,size(dat,2)]);
    % z-score units
    % use averaged std to account for zero spike baselines!
    base = std(dat(:,blank),[],2);
    base = ones(size(base))*nanmean(base);
    if sum(base == 0)
      % then units will fail
      base = std(dat(:,blank),[],2);
      base = ones(size(base))*prctile(base,75);
    end
   if sum(base == 0 )
     base = ones(size(base));
   end
   dat = dat./repmat(base,[1,size(dat,2)]);
    
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functiont time-adjust individual trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  tmp = subDataOffsetCorrect(tmp,offset)

d = size(tmp,2);
for k=1:size(tmp,1)
  y = tmp(k,:);
  y([1:d-offset(k)]) = tmp(k,[offset(k)+1:end]);
  tmp(k,:) = y;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to smooth single trial spike trains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out =  sub_smooth_spk(Matrix,width)

if rem(width,2)
  width = width+1;
end
gauss = gausswin(width);
gauss = gauss / sum(gauss);
for k=1:size(Matrix,1)
  resp = Matrix(k,:);
  resp2 = conv(resp,gauss);
  resp2 = resp2( floor(length(gauss)/2):end);
  resp2 = resp2(1:length(resp));
  out(k,:) = resp2;
end
return;




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specialized function for the adaptation paradigm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Resp = SubGetAdapt(Data,ARG)
WIN_PRE = 150;
WIN_POST = 500;
    
% do a trivial time adjustment
ARG.etc.time_adjust = 0; 
[time_adjust,computed_adjust] = subGetTrialTiming(Data,ARG);
adjust.time_adjust = time_adjust;
adjust.computed_adjust = computed_adjust;
Resp.adjust = adjust;

% compute variability of signal to remove (bad) trials
if ARG.etc.bad_trial_thr_lfp
  Resp.data_quality = subGetDataQ(Data,ARG,time_adjust);
else
  Resp.data_quality.mask =[];
  Resp.data_quality.num_bad = zeros(length(Data),length(Data{1}.info.BEH.conds));
end
MASK = Resp.data_quality.mask;
if isempty(MASK)
  do_mask = 0;
else
  do_mask = 1;
end
% ....................................................................
% obtain a list of conditions and stimuli
X = Data{1}.info.BEH.RevCor;
tmp =[];
for k=1:length(X), tmp = [tmp; X{k}]; end;
conds = unique(tmp(:,2));

% so far we do only the Amua and Spk for this paradigm
% ......................................................................
% AMUA
if ~isfield(Data{1},'Amua')
  fprintf('Amua not found'); 
else
  fprintf('Amua\t');
  Nel = length(Data);
  rate = Data{1}.info.ARG.Amua.rate;
  for E=1:Nel
    clear Timing dat
    dat = Data{E}.Amua; 
    % add time-adjust to individual trials
    offset = round(time_adjust{1}/1000*rate);
    dat = dat{1}; % we assume that there is only one stimulus
    dat = subDataOffsetCorrect(dat,offset);
    Nstim = length(dat);
    t_ax = Data{E}.Amua_taxis;
    T_B = find(t_ax<=-50);
    % remove bad trials
    if do_mask
      good = find(MASK{E,1}==0);
      dat = dat(good,:);
      for c=1:length(good)
        Timing{c} = Data{E}.info.BEH.RevCor{good(c)};
      end
    else
      Timing = Data{E}.info.BEH.RevCor;
    end    
    dat = subDataPrep(dat,ARG,T_B,'Amua');  
    % separate the different categories:
    Nstim = size(Timing{1},1);
    for C=1:length(conds)
      CondD{C} = [];
      for S=1:Nstim
        CondReps{C,S} =[];
      end
    end
    WINL = length(find( (t_ax>=0-WIN_PRE).*(t_ax<=0+WIN_POST)));
    for T=1:length(Timing)
      C = Timing{T}(2,2)+1;
      % store evp
      CondD{C} = [CondD{C};dat(T,:)];
      % compute amplitude for each stimulus (here the timing across trials
      % is more precise
      for S=1:Nstim
        [dummy,ton] = min(abs((Timing{T}(S,4)-Timing{T}(1,4)-WIN_PRE)-t_ax));
        T_INT = ton+[1:WINL];
        CondReps{C,S} = [CondReps{C,S}; dat(T,T_INT)];
      end
    end
    for C=1:length(conds)
      Amua{E}.Evp{C} = median(CondD{C},1);
      % test for responsiveness in each condition
      evp2 = subDataPrep(mean(CondD{C},1),ARG,T_B,'Amua');
      Amua{E}.Evp_sigrespZ(C,:,:) = ck_physio_getlatency(evp2,t_ax);    
      Amua{E}.t_axis = Data{E}.Amua_taxis;
      % each stimulus
      for S=1:size(CondReps,2)
        Amua{E}.Resp{C,S} = CondReps{C,S};
        Amua{E}.Resp_tax = t_ax(find( (t_ax>=0-WIN_PRE).*(t_ax<=0+WIN_POST)));
      end
    end
    Amua{E}.elec = Data{E}.elec;
    Amua{E}.grpname = Data{E}.info.grpname;
    Amua{E}.site = Data{E}.info.site;
  end % E
  Resp.Amua = Amua;
end
% ......................................................................
% SPK
if ~isfield(Data{1},'Spk')
  fprintf('Spks not found'); 
  for k=1:length(Resp.Amua)
    Resp.Spk{k}.nunits = 0;
  end
else
  fprintf('Spk\t');
  Nel = length(Data);
  bin = ARG.Spk.bin;
  for E=1:Nel
    clear Timing dat
    dat = Data{E}.Spk;
    if ~isempty(dat)
      trial_dur = ceil( abs(diff(Data{E}.Lfp_taxis([1,end]))) / bin)+5;
      Nunits = length(dat{1});
      for U=1:Nunits
        tmp = dat{1}{U};
        Ntrial = length(tmp);
        Matrix = zeros(Ntrial,trial_dur);
        for t=1:Ntrial
          if ~isempty(tmp{t})
            spk_time = tmp{t}(:,1)*1000; % units of ms
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % time-adjust spike timing !!!!!!!!!!!!!!!!
            % this  might create spikes before trial onset, remove these
            spk_time = spk_time - time_adjust{1}(t);  
            spk_time = ceil(spk_time / bin); % binning
            spk_time = spk_time(find(spk_time>0));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for l=1:length(spk_time)
              Matrix(t,spk_time(l)) = Matrix(t,spk_time(l)) + 1;
            end
          end % if
        end % trial
        Matrix = Matrix * (1000/bin); % THIS ONLY FOR UNITS OF SPK/sec
        % remove bad trials
        if do_mask
          good = find(MASK{E,1}==0);
          Matrix = Matrix(good,:);
          for c=1:length(good)
            Timing{c} = Data{E}.info.BEH.RevCor{good(c)};
          end
        else
          Timing = Data{E}.info.BEH.RevCor;
        end        
        Nstim = size(Timing{1},1);
        for C=1:length(conds)
          CondD{C} = [];
          for S=1:Nstim
            CondReps{C,S} =[];
          end
        end
        % smooth
        width = ARG.Spk.psth_smooth/ARG.Spk.bin;
        width = ceil(width*1.8);
        Matrix = sub_smooth_spk(Matrix,width);
        model = zeros(1,size(Matrix,2));
        t_ax = [1:length(model)] -1 + round(Data{E}.Spk_taxisoffset/bin);
        t_ax = t_ax*bin;
        
        WINL = length(find( (t_ax>=0-WIN_PRE).*(t_ax<=0+WIN_POST)));
        for T=1:length(Timing)
          C = Timing{T}(2,2)+1;
          % store evp
          CondD{C} = [CondD{C}; Matrix(T,:)];
          % compute amplitude for each stimulus (here the timing across trials
          % is more precise
          for S=1:Nstim
            [dummy,ton] = min(abs((Timing{T}(S,4)-Timing{T}(1,4)-WIN_PRE)-t_ax));
            T_INT = ton+[1:WINL];
            CondReps{C,S} = [CondReps{C,S}; Matrix(T,T_INT)];
          end
        end % T
        
        for C=1:length(conds)
          SPK{E}.Mean_resp{U}{C} = mean(CondD{C},1);
          SPK{E}.Mean_resp_tax = t_ax;
          for S=1:size(CondReps,2)
            SPK{E}.Resp{U}{C,S} = CondReps{C,S};
            SPK{E}.Resp_tax = t_ax(find( (t_ax>=0-WIN_PRE).*(t_ax<=0+WIN_POST)));
          end
        end % C
        SPK{E}.Obs_id{U} = single(Data{E}.Lfp_obsid{1}(good,:));
        Wave = Data{E}.Spk_detail{U}(:,8:end-1);
        if ~isempty(Wave)
          xmean = mean(Wave,1);
          sem = std(Wave)/sqrt(size(Wave,1));
          crit = tinv((1 - 0.05 / 2), (size(Wave,1)-1)) .* sem;
          SPK{E}.wave{U}(1,:) = xmean;
          SPK{E}.wave{U}(2,:) = xmean - crit;
          SPK{E}.wave{U}(3,:) = xmean + crit;
          SPK{E}.wave{U}(4,:) = std(Wave);
          % signal to noise ratio of wave form
          b = std(xmean([1:5,21:22]));
          p = max(abs(xmean));
          SPK{E}.wave_snr(U) = p/b;
          % ISI histogram with 0.5 ms bin
          isi = Data{E}.Spk_detail{U}(:,3);
          isi = diff(isi*1000);
          isi = hist(isi,[0:0.5:400]);
          SPK{E}.wave_ISI{U} = isi/sum(isi);
        else
          SPK{E}.wave{U} = [];
          SPK{E}.wave_snr(U) = 0;
          SPK{E}.wave_ISI{U} =[];
        end
        
      end % unit
      SPK{E}.nunits =  Nunits;
    else
      SPK{E}.nunits = 0;
    end % if unit
  end % E
  Resp.SPK = SPK;
end % if spike



return;


