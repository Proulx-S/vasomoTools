function J = compMTall(tapers,data,Fs,tvec)

data=change_row_to_column(data);
[NC,C]=size(data); % size of data
[NK,K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
NT = NC; clear NC NK
% % tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])
% tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
% % tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])
% data=data(:,:,ones(1,K)); % add taper indices to data
% % tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])
% data=permute(data,[1 3 2]); % [time X tapers X vox] [NC K C] reshape data to get dimensions to match those of tapers
% % tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])
% data=data.*tapers; clear tapers % [time X tapers X vox] [NC K C] product of data with tapers
% % tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])

%% Project data in taper space
tic
tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])
data=permute(data(:,:,ones(1,K)),[1 3 2]).*tapers(:,:,ones(1,C)); clear tapers
toc
tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])

varList = {{'data'}}; for i = 1:length(varList); eval(['tmp = whos(''' strjoin(varList{i},''',''') ''');']); eval(['Gb.' strjoin(varList{i},'') ' = sum([tmp.bytes])./1e9;']); end


%% Compute coherence
if ~exist('tvec','var') || isempty(tvec)
    tvec=1/Fs *(0:NT-1)';
end
tvec=tvec*2*pi*1i;
nfft=max(2^(nextpow2(NT)+0),NT);
freqs=getfgrid(Fs,nfft,[0 Fs/2]);
NF=length(freqs);
s = nan(1,K,NF,'single');

data
tapers
tvec
data = permute(data,[2 1 3 4]);
tapers = permute(tapers,[1 2 3 4]);
freqs = permute(freqs,[1 3 2 4]);
tvec = permute(tvec,[4 1 2 3]);
whos data tapers tvec freqs

nf = NF;
fInd = 1:nf;
tic
[~,s,~] = pagesvd(reshape(data*reshape(tapers.*exp(-freqs(fInd).*tvec),[NT K*nf]),[C K nf]),"econ","vector");
toc
whos u s v


% tfvec = exp(-freqs.*tvec);
% whos data tapers tfvec
% size(data*reshape(tapers.*tfvec,[NT K*NF]))
% size(reshape(data*reshape(tapers.*tfvec,[NT K*NF]),[C K NF]))

% s = nan(1,K,NF,'single');
% nf = 1;
% tic
% for fInd = 1:10
%     % size(reshape(data*reshape(tapers.*exp(-freqs(fInd).*tvec),[NT K*nf]),[C K nf]))
%     [u,s(:,:,fInd),v] = compSVD(reshape(data*reshape(tapers.*exp(-freqs(fInd).*tvec),[NT K*nf]),[C K nf]));
% end
% toc


% [u,s(:,:,fInd),v] = compSVD(reshape(data*reshape(tapers.*exp(-freqs(fInd).*tvec),[NT K*nf]),[C K nf]));
toc


freqs = permute(-freqs,[4 2 1 3]);
tvec = permute(tvec,[4 1 2 3]);
tapers = permute(tapers,[1 3 2]);
whos freqs tvec tapers

projBloc = permute(freqs(1:10),[1 2 3 4])...
         .*...
         permute(tvec,[1 2 3 4])...
         .*...
         permute(tapers,[1 2 3 4]);

projBloc

fvec
exp(-freqs(fInd).*squeeze(tvec))

NFb = 10;
fInd = (1:NFb)+1;
tic
tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])
[u,s(:,:,fInd),v] = compSVD(permute(sum(data.*squeeze(exp(-freqs(fInd).*squeeze(tvec))),1),[3 2 1]));
tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])
toc


parfor fInd = 1:NF
    [u,s(:,:,fInd),v] = compSVD(permute(sum(data.*squeeze(exp(-freqs(fInd)*tvec)),1),[3 2 1]));
end
c = permute(s/sum(s,2),[2 1]); % [vox X tapers X freq] [C K NF]


varList = {'u' 's' 'v' 'c'}; eval(['tmp = whos(''' strjoin(varList,''',''') ''');']); eval(['Gb.' strjoin(varList,'') ' = sum([tmp.bytes])./1e9;']);
disp(['cost of ' strjoin(varList,'') ' is ' num2str(Gb.(strjoin(varList,''))./Gb.data) ' the size of data'])
toc
tmp = whos; disp(['total mem: ' num2str(sum([tmp.bytes])/1e9) 'GB'])

% % !!!! Grrr, need the all tapers at a given frequency
% tvec=1/Fs *(1:N)';
% % tvec=repmat(tvec,[1 K]);
% tvec=tvec*2*pi*1i;
% freqs=getfgrid(Fs,nfft,fpass);
% nf=length(freqs);
% if nf>1; error('something wrong'); end
% if ~isempty(mdkp) && mdkp > min(K,NCHAN); error('mdkp has to be less than both K and NCHAN'); end
% tapersXfreqs=tapers.*exp(-freqs*tvec);


%% Compute fft
tic
nfft=max(2^(nextpow2(NT)+pad),size(data,1));
%%% in blocs
J=fft(data,nfft)/Fs; %clear data   % fft of projected data
varList = {'J'}; eval(['tmp = whos(''' strjoin(varList,''',''') ''');']); eval(['Gb.' strjoin(varList,'') ' = sum([tmp.bytes])./1e9;']);
disp(['cost of ' strjoin(varList,'') ' is ' num2str(Gb.(strjoin(varList,''))./Gb.data) ' the size of data'])
toc

%% Compress fft results
[f,findx]=getfgrid2(Fs,nfft,fpass); 
J = fftshift(J,1);
J=J(findx,:,:);
J = permute(J,[1 3 2]);

%% Average fft powers across tapers
J = sqrt(mean(abs(J).^2,2));

function [u,s,v] = compSVD(data)
[u,s,v] = svd(data,"econ"); % [vox X tapers] [C K]
s = diag(s).^2;

 





% dist = repmat(round(C/nBloc),[1 nBloc]); dist(end) = dist(end) + (C - sum(dist));
% data = mat2cell(data,NC,dist);
% J = cell(size(data));
% for bloc = 1:nBloc
%     tic
%     disp([num2str(bloc) '/' num2str(nBloc) ': computing'])
%     C = dist(bloc);
%     J{bloc} = computeVectorized(tapers,data{bloc},NC,K,C,nfft,Fs);
%     % average power to save space (phase information is lost)
%     J{bloc} = sqrt(mean(conj(J{bloc}).*J{bloc},2));
%     % delete data to save space
%     data{bloc} = {};
%     disp([num2str(bloc) '/' num2str(nBloc) ': done'])
%     toc
% end
% J = cat(3,J{:});
