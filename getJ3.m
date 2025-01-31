function [JChronux,JFFT] = getJ3(d,tp,t,f,noFFT)
% [time x trial x run x taper x freq x vox x window]
%[N ,E ,R ,K,F,V,W]
if ~exist('noFFT','var') noFFT = []; end
if isempty(noFFT);       noFFT = 0; end


 [Nd,Ed,Rd,~,~,V,Wd] = size(d);
 [~ ,~ ,~ ,K,~,~,~ ] = size(tp);
 [~ ,~ ,~ ,~,F,~,~ ] = size(f);
 [Nt,Et,Rt,~,~,~,Wt] = size(t);
if Nd~=Nt         ; dbstack; end; N = Nd;
if Et~=1 && Ed~=Et; dbstack; end; E = Ed;
if Rt~=1 && Rd~=Rt; dbstack; end; R = Rd;
if Wt~=1 && Wd~=Wt; dbstack; end; W = Wd;

JChronux = zeros(1,E,R,K,F,V,W);
JFFT     = zeros(1,E,R,K,F,V,W);
% [N ,E ,R ,K,F,V,W]
for e = 1:E
    % try
    if size(tp,2)==E
        tp2 = reshape(  tp(:,e,:,:,:,:,:,:) .* exp(-f.*t(:,e)*2*pi*1i)                   ,[N 1*R*K*F*1*Wt]);
    else
        tp2 = reshape(  tp(:,1,:,:,:,:,:,:) .* exp(-f.*t(:,e)*2*pi*1i)                   ,[N 1*R*K*F*1*Wt]);
    end
        d2  = reshape(  d(:,e,:,:,:,:,:,:)                              ,[N 1*R*1*1*V*W ]);
        % d2  = reshape(  d(:,e,:,:,:,:,:,:) - mean(d(:,e,:,:,:,:,:,:),1) ,[N 1*R*1*1*V*W ]); % removing the mean here

        % tp2 = reshape(  tp .* exp(-f.*t*2*pi*1i)  ,[N Et*R*K*F*1*Wt]);
        % d2  = reshape(  d - mean(d,1)             ,[N E*R*1*1*V*W]); % removing the mean here

        d2 = permute(d2,[2 1]);
        j = d2*tp2; % [V 1*R*K*F*1*W]
        % [E*R*1*1*V*W Et*R*K*F*1*Wt]

        j = reshape(j,[1*R*1*1*V*W 1 R K F 1 Wt]);
        j = permute(j,[2 3 4 5 6 7 1]);
        j = reshape(j,[1 R K F 1 Wt 1 R 1 1 V W]);
        % jChronux = permute(j,[13 7 2 3 4 11 6 1 5 8 9 10 12]);
        JChronux(:,e,:,:,:,:,:) = permute(j,[13 7 2 3 4 11 6 1 5 8 9 10 12]); %[N ,E ,R ,K,F,V,W]
    %     if any(size(J,9:15)~=1); dbstack; error('X'); end
    % 
    %     if any([E R W]>1)
    %         dbstack;
    %         error('double-check that')
    %     end
    % 
    % catch err
    %     % try a low memory approach
    %     % d2 = reshape(  d(:,e,:,:,:,:,:,:)                              ,[N 1*R*1*1*V*W ]);
    %     % d2 = permute(d2,[2 1]);
    % 
        if any([R W]>1)
            dbstack;
            error('double-check that')
        end
        NFFT = (F-1)*2;
        d2 = reshape(d(:,e,:,:,:,:,:,:),[N 1*R*1*1*V*W ]); % [N ,E ,R ,K,F,V,W] -> [N,E*R*K*F*V*W]
        j = fft(d2.*tp,NFFT,1); % [N,E*R*K*F*V*W].*[N,E,R,K] -> [NFFT,Ed*Rd*Kd*Fd*Vd*Wd,Rtp,Ktp]
        j = j(1:NFFT/2+1,:,:,:,:,:,:,:);% [Ff,Ed*Rd*Kd*Fd*Vd*Wd,Rtp,Ktp]      [F V 1 K]
        j = permute(j,[1 3 4 2]); % [Ff,Ed*Rd*Kd*Fd*Vd*Wd,Rtp,Ktp] -> [Ff,Rtp,Ktp,Ed*Rd*Kd*Fd*Vd*Wd];
        j = reshape(j,[F 1 K 1 R 1 1 V W]); %[Ff,Rtp,Ktp,Ed,Rd,Kd,Fd,Vd,Wd]
        % jFFT = permute(j,[2 4 5 3 1 8 9 6 7]);
        JFFT(:,e,:,:,:,:,:) = permute(j,[2 4 5 3 1 8 9 6 7]); %[F,Rtp,K,Ed,R,Kd,Fd,V,W]
        
        % % % j = permute(j,[2 1 3 5 4 6 7 8]);
        % % % jFFT = permute(j,[3 6 4 5 2 1 8 7]);
        % % % J(:,e,:,:,:,:,:) = permute(j,[3 6 4 5 2 1 8 7]); %[N ,E ,R ,K,F,V,W]
        % % size(Y)
        % % load tmp j
        % % size(j)
        % [N ,E ,R ,K,F,V,W]
        % size(d)
        % {size(jChronux) size(jFFT)}'
        % size(f)
        % 
        % % min(abs(jChronux(1,1,1,1,:,1)) - abs(jFFT(1,1,1,1,:,1)))
        % % max(abs(jChronux(1,1,1,1,:,1)) - abs(jFFT(1,1,1,1,:,1)))
        % max(abs(abs(  jChronux(1,1,1,1,:,1)) - abs(  jFFT(1,1,1,1,:,1))))
        % max(abs(angle(jChronux(1,1,1,1,:,1)) - angle(jFFT(1,1,1,1,:,1))))
        % 
        % figure('WindowStyle','docked');
        % plot(squeeze(f),squeeze(real(jChronux(1,1,1,1,:,1))),'r'); hold on
        % plot(squeeze(f),squeeze(real(jFFT(1,1,1,1,:,1))),'--k')
        % plot(squeeze(f),squeeze(imag(jChronux(1,1,1,1,:,1))),'g'); hold on
        % plot(squeeze(f),squeeze(imag(jFFT(1,1,1,1,:,1))),'--k')
        % % 
        % % 
        % % 
        % % 
        % % 
        % squeeze(jChronux(1,1,1,1,:,1))
        % squeeze(jFFT(1,1,1,1,:,1))'
        % 
        % figure('WindowStyle','docked');
        % x = abs(squeeze(jChronux(1,1,1,1,:,1)));
        % % x = x./mean(x);
        % y = abs(squeeze(jFFT(1,1,1,1,:,1)))';
        % max(abs(x-y))
        % % y = y./mean(y);
        % % imagesc((x - y)./(x + y)); colorbar
        % imagesc(abs(x - y)); colorbar
        % ax = gca; ax.DataAspectRatio = [1 1 1];
        % ax.ColorScale = 'log';
        % 
        % 
        % k = 1;
        % figure('WindowStyle','docked');
        % plot(squeeze(f),squeeze(abs(jChronux(1,1,1,k,:,1))),'k'); hold on
        % plot(squeeze(f),squeeze(abs(jFFT(1,1,1,k,:,1))),'--k','LineWidth',5)
        % yyaxis right
        % plot(squeeze(f),squeeze(angle(jChronux(1,1,1,k,:,1))),'r'); hold on
        % plot(squeeze(f),squeeze(angle(jFFT(1,1,1,k,:,1))),'--r','LineWidth',5)
        % ax = gca; ax.YAxis(2).Color = 'r';
        % xlim([0 0.1])
    % 
    % end

    
end
if ~noFFT
    % keyboard
end
% J = JChronux;
% J = reshape(J,[V E R K F 1 W]); % [V E R K F 1 W]
% J = permute(J,[6 2 3 4 5 1 7]); % [N E R K F V W]