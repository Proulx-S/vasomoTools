function J = getJ4(d,tp,t,f,chronuxFlag,testFlag)
% [time x trial x run x taper x freq x vox x window]
%[N ,E ,R ,K,F,V,W]
if ~exist('chronuxFlag','var'); chronuxFlag = []; end
if ~exist('testFlag','var');       testFlag = []; end
if isempty(chronuxFlag);        chronuxFlag = 0; end
if isempty(testFlag);              testFlag = 0; end


[Nd,Ed,Rd,~,~,V,Wd] = size(d);
[~ ,~ ,~ ,K,~,~,~ ] = size(tp);
[~ ,~ ,~ ,~,F,~,~ ] = size(f);
[Nt,Et,Rt,~,~,~,Wt] = size(t);
if Nd~=Nt         ; dbstack; end; N = Nd;
if Et~=1 && Ed~=Et; dbstack; end; E = Ed;
if Rt~=1 && Rd~=Rt; dbstack; end; R = Rd;
if Wt~=1 && Wd~=Wt; dbstack; end; W = Wd;

if chronuxFlag
    J = doItChronuxStyle(d,tp,t,f,N,E,R,K,F,V,W,Wt);
else
    NFFTexpected = 2.^(nextpow2(N) + (0:10));
    if ~ismember(2*(length(f)-1),NFFTexpected)
        dbstack; error(['using fftStyle computation of J, but requested frequencies do not satisfy fft requirements' newline 'turn on chronuxFlag for getJ4.m']);
    end
    J = doItFFTStyle(d,tp,t,f,N,E,R,K,F,V,W,Wt);
end

if testFlag && ~chronuxFlag
    try
        tic
        JC = doItChronuxStyle(d,tp,t,f,N,E,R,K,F,V,W,Wt);
        toc

        {size(JC) size(J)}'

        max(abs(           abs(JC(:))   - abs(J(:))) )
        max(abs(wrapToPi(angle(JC(:)) - angle(J(:)))))

        % figure('WindowStyle','docked');
        % plot(squeeze(f),squeeze(abs(JC(1,1,1,1,:,1))))
        % hold on
        % % yyaxis right
        % plot(squeeze(f),squeeze(abs(J(1,1,1,1,:,1))))
        % legend({'JC' 'JF'})
        %
        % figure('WindowStyle','docked');
        % plot(squeeze(f),squeeze(angle(JC(1,1,1,1,:,1))))
        % hold on
        % plot(squeeze(f),squeeze(angle(J(1,1,1,1,:,1))))
        % legend({'JC' 'JF'})
    catch ME
        JC = [];
        warning([ME.message newline 'skipping this'])
    end

end


function J = doItChronuxStyle(d,tp,t,f,N,E,R,K,F,V,W,Wt)
J = zeros(1,E,R,K,F,V,W);
% [N ,E ,R ,K,F,V,W]
for e = 1:E
    if size(tp,2)==E
        tp2 = reshape(  tp(:,e,:,:,:,:,:,:) .* exp(-f.*t(:,e,:,:,:,:,:,:)*2*pi*1i) ,[N 1*1*K*F*1*Wt*1]); % [N,E*R*K*F*V*W*M]
    else
        tp2 = reshape(  tp(:,1,:,:,:,:,:,:) .* exp(-f.*t(:,e,:,:,:,:,:,:)*2*pi*1i) ,[N 1*1*K*F*1*Wt*1]); % [N,E*R*K*F*V*W*M]
    end
    d2  = reshape(   d(:,e,:,:,:,:,:,:)                                        ,[N 1*R*1*1*V*W *1 ]);% [N,E*R*K*F*V*W*M]
    d2 = permute(d2,[2 1]);
    % [E*R*K*F*V*W*M,N]
    % [1*R*1*1*V*W*1,N]
    tp2;
    % [N,E*R*K*F*V*W *M]
    % [N,1*1*K*F*1*Wt*1]
    j = d2*tp2;
    % [E*R*K*F*V*W*M,E*R*K*F*V*W*M]
    % [1*R*1*1*V*W*1,1*1*K*F*1*W*1]
    j = reshape(j,[1*R*1*1*V*W*1 1 1 K F 1 W 1]);
    j = permute(j,[2 3 4 5 6 7 8 1]);
    j = reshape(j,[1 1 K F 1 W 1 1 R 1 1 V W 1]);
    % [E R K F V W M E R K F V W M]
    % [1 1 K F 1 W 1 1 R 1 1 V W 1]
    % [tp2           d2           ]
    j = permute(j,[1 8 9 3 4 12 13 7 2 5 6 10 11 14 15 16]);
    % [N,E,R,K,F,V,W,M]
    % [1,1,R,K,F,V,W,1]

    J(:,e,:,:,:,:,:) = j;
end



function J = doItFFTStyle(d,tp,t,f,N,E,R,K,F,V,W,Wt)
J = zeros(1,E,R,K,F,V,W);

for e = 1:E
    if any([R W]>1)
        dbstack;
        error('double-check that')
    end
    NFFT = (F-1)*2;
    d2 = reshape(d(:,e,:,:,:,:,:,:),[N 1*R*1*1*V*W ]); % [N ,E ,R ,K,F,V,W] -> [N,E*R*K*F*V*W]
    if size(tp,2)==E
        j = fft(d2.*tp(:,e,:,:,:,:,:,:),NFFT,1); % [N,E*R*K*F*V*W].*[N,E,R,K] -> [NFFT,Ed*Rd*Kd*Fd*Vd*Wd,Rtp,Ktp]
    else
        j = fft(d2.*tp(:,1,:,:,:,:,:,:),NFFT,1); % [N,E*R*K*F*V*W].*[N,E,R,K] -> [NFFT,Ed*Rd*Kd*Fd*Vd*Wd,Rtp,Ktp]
    end
    j = j(1:NFFT/2+1,:,:,:,:,:,:,:);% [Ff,Ed*Rd*Kd*Fd*Vd*Wd,Rtp,Ktp]      [F V 1 K]
    j = permute(j,[1 3 4 2]); % [Ff,Ed*Rd*Kd*Fd*Vd*Wd,Rtp,Ktp] -> [Ff,Rtp,Ktp,Ed*Rd*Kd*Fd*Vd*Wd];
    j = reshape(j,[F 1 K 1 R 1 1 V W]); %[Ff,Rtp,Ktp,Ed,Rd,Kd,Fd,Vd,Wd]
    J(:,e,:,:,:,:,:) = permute(j,[2 4 5 3 1 8 9 6 7]); %[F,Rtp,K,Ed,R,Kd,Fd,V,W]

    % rotate to phase=0 at t=0
    J(:,e,:,:,:,:,:) = J(:,e,:,:,:,:,:).*exp(-f.*t(1,e,:,:,:,:,:,:)*2*pi*1i);
end
