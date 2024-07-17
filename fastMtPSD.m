function [J,f] = fastMtPSD(tapers,data,Fs,Fpass,nBloc,avTapers,verbose,tvec)
% tapers and tvec must be [time x trial x tapers     ]
% data must be            [time x vox   x trial x run]
if ~exist('tvec','var');   tvec = []; end
if ~exist('Fpass','var'); Fpass = []; end
if isempty(Fpass);        Fpass = [0 inf]; end

[N,C,E,R]  = size(data);
[Nk,Ek,K]  = size(tapers);
if ~isempty(tvec)
    [Nv,Ev,Kv] = size(tvec);
    if ~all([N*E==Nk E==Ek Nk==Nv Ek==Ev K==Kv]); dbstack; error('inappropriate size for data, tapers or tvec'); end
else
    if ~all([N*E==Nk E==Ek]); dbstack; error('inappropriate size for data, tapers or tvec'); end
end
% NE = N*E;
% sz     = size(data,1:4);
% NR     = sz(4);
% NE = sz(3);
% NC     = sz(1)*NE;
% C = prod(sz([2 4]));
% [NC,C]=size(data); % size of data
% if NE~=Nk; dbstack; error('length of tapers is incompatible with length of data'); end
% if NK~=NC; error('length of tapers is incompatible with length of data'); end
% NT = NC; clear NC NK

pad = 0;
NFFT=max(2^(nextpow2(N*E)+pad),N*E);
% NFFT=max(2^(nextpow2(NT)+pad),NT);
if Fpass(2)==inf; Fpass(2) = Fs/2; end
[f,fInd]=getfgrid(Fs,NFFT,Fpass);
F = length(f);

if ~exist('nBloc','var') || isempty(nBloc)
    nBloc = 1;
end
if ~exist('avTapers','var') || isempty(avTapers)
    avTapers = false;
    if nBloc>1
        avTapers = true;
    end
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = 2;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute everything at once for speed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nBloc==1
    if isempty(tvec)
        %% Using FFT
        switch avTapers
            case 0
                % Project to taper space -> compute fft -> keep data from all tapers
                J =            reshape(  fft(reshape( data(:,:).*tapers ,[N C*K]),NFFT)/Fs  ,[NFFT C K])         ;
            case 1
                if verbose>1; disp('averaging across tapers--phase information will be lost'); end
                % Project to taper space -> compute fft -> convert to power -> average across tapers
                J = mean(abs(  reshape(  fft(reshape( data(:,:).*tapers ,[N C*K]),NFFT)/Fs  ,[NFFT C K])  ).^2,3);
            case 2
                if verbose>1; disp('phase-coherent averaging across tapers'); end
                % Project to taper space -> compute fft
                J = reshape(  fft(reshape( data(:,:).*tapers ,[N C*K]),NFFT)/Fs  ,[NFFT C K]);
                %%%%%%%%%
                % align phases across trials
                dbstack; error('need to align phases across trials here')
                %%%%%%%%%
                % average across tapers
                J = abs(mean(  J  ,3)).^2;
            otherwise
        end
        % Keep one-sided spectrum
        J = permute(J(fInd,:,:),[1 3 2]); % [freq x taper x vox]

    else
        %%% Use fourrier series--this allows manipulations of tvec to align phase across trials
        if R>1; dbstack; error('double-check that'); end
        
        tvec2 = exp(  permute(f,[1 3 2])  .*  reshape(tvec*2*pi*1i,[N*E E*K])  ); %[N*E x E*K   x F]
        tapers2 = reshape(  reshape(tapers,[N*E E*K])  .*  tvec2  ,[N*E E*K*F]);  %[N*E x E*K*F    ]
        data2 = reshape(permute(data,[2 1 3]),[C N*E]);
        J = permute(  reshape(  data2 * tapers2  ,[C E K F])  ,[4 3 1 2]) / Fs; % [freq x taper x vox x trial]

        switch avTapers
            case 0 % don't average across tapers (or trials)
            case 1 % average across tapers (and trials)
                J = mean(  mean( abs(J).^2 ,4)  ,2);
            case 2 % phase-coherent average across trials and average across tapers
                J = mean(  abs( mean(J,4) ).^2  ,2);
        end
        
        % switch avTapers
        %     case 0 % don't average across tapers
        %         J =             permute(  reshape(  permute(data,[2 1]) * reshape(tapers.*exp(-f.*permute(tvec*2*pi*1i,[1 3 2])),[NT F*K])  ,[C F K])/Fs  ,[2 3 1])          ;
        %     case 1 % average across tapers
        %         J = mean(abs(   permute(  reshape(  permute(data,[2 1]) * reshape(tapers.*exp(-f.*permute(tvec*2*pi*1i,[1 3 2])),[NT F*K])  ,[C F K])/Fs  ,[2 3 1])   ).^2,2);
        %     case 2 % phase-coherent average across tapers
        %         J = abs(mean(   permute(  reshape(  permute(data,[2 1]) * reshape(tapers.*exp(-f.*permute(tvec*2*pi*1i,[1 3 2])),[NT F*K])  ,[C F K])/Fs  ,[2 3 1])   ,2)).^2;
        % end
        % hold off
        % for kInd = 1:3
        %     fInd = 6;
        %     t = tvec(1,kInd):1/Fs:tvec(end,kInd);
        %     fSeries = nan(size(t,2),1);
        %     fSeries(ismembertol(t,tvec(:,kInd)')) = tapers(:,1,kInd).*exp(-f(fInd).*permute(tvec(:,kInd)*2*pi*1i,[1 3 2]));
        %     plot(t,real(fSeries)); ylim([-0.5 0.5]); hold on
        % end


        % if avTapers
        %     J = mean(abs(   permute(  reshape(  permute(data,[2 1]) * reshape(tapers.*exp(-f.*squeeze(tvec*2*pi*1i)),[NT NF*K])  ,[C NF K])/Fs  ,[2 3 1])   ).^2,2);
        % else
        %     J =             permute(  reshape(  permute(data,[2 1]) * reshape(tapers.*exp(-f.*squeeze(tvec*2*pi*1i)),[NT NF*K])  ,[C NF K])/Fs  ,[2 3 1])          ;
        % end
    end
else
    if ~isempty(tvec); dbstack; error('code that'); end
    % Split computations across nBloc blocs of roughly equal number of
    % voxels, averaging across tapers and deleting data between each bloc
    % to save on memory
    if verbose>1; disp('averaging across tapers--phase information will be lost'); end
    Cdist = repmat(round(C/nBloc),[1 nBloc]); Cdist(end) = Cdist(end) + (C - sum(Cdist));
    data = mat2cell(data(:,:),NT,Cdist);
    J = cell(size(data));
    for bloc = 1:nBloc
        tic
        if verbose; disp(['bloc ' num2str(bloc) ' of ' num2str(nBloc) ': computing']); end
        % Project to taper space -> compute fft -> convert to power ->
        % average across tapers -> sqrt for consistency.
        J{bloc}=sqrt(mean(abs(reshape(fft(reshape(data{bloc}.*tapers,[NT Cdist(bloc)*K]),NFFT)/Fs,[NFFT Cdist(bloc) K])).^2,3));
        toc

        % Delete data to save space
        data{bloc} = {};
        if verbose; disp(['bloc ' num2str(bloc) ' of ' num2str(nBloc) ': done']); end
        toc
    end
    % Catenate compute blocs
    J = cat(2,J{:});
end


% %% Adjust spectrum and keep data with Fpass
% df=Fs/NFFT;
% f=0:df:(Fs-df);
% f = f-Fs/2;
% f = -f;
% f = fftshift(f);
% % fInd = find(f>=Fpass(1) & f<=Fpass(2));
% % [~,b] = sort(f(fInd));
% % fInd = fInd(b);
% if Fpass(2)==inf; Fpass(2) = Fs/2; end
% 
% [f,fInd]=getfgrid(Fs,NFFT,Fpass);
% 
% sz = [length(fInd) K C/NR NR];
% if avTapers; sz(2) = 1; end
% 
% if isempty(tvec)
%     J = reshape(permute(J(fInd,:,:),[1 3 2]),sz); % [freq x 1 x vox]
% else
%     J = reshape(permute(J,[1 3 2]),sz);
% end
% % f = f(fInd);