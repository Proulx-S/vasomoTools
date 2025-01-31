function J = getJ2(d,tp,t,f)
% [time x trial x run x taper x freq x vox x window]
%[N ,E ,R ,K,F,V,W]
 [Nd,Ed,Rd,~,~,V,Wd] = size(d);
 [~ ,~ ,~ ,K,~,~,~ ] = size(tp);
 [~ ,~ ,~ ,~,F,~,~ ] = size(f);
 [Nt,Et,Rt,~,~,~,Wt] = size(t);
if Nd~=Nt         ; dbstack; end; N = Nd;
if Et~=1 && Ed~=Et; dbstack; end; E = Ed;
if Rt~=1 && Rd~=Rt; dbstack; end; R = Rd;
if Wt~=1 && Wd~=Wt; dbstack; end; W = Wd;

% for e = 1:E
% tp2 = reshape(  tp .* exp(-f.*t(:,e)*2*pi*1i)                   ,[N 1*R*K*F*1*Wt]);
% d2  = reshape(  d(:,e,:,:,:,:,:,:) - mean(d(:,e,:,:,:,:,:,:),1) ,[N 1*R*1*1*V*W ]); % removing the mean here
tp2 = reshape(  tp .* exp(-f.*t*2*pi*1i)  ,[N Et*R*K*F*1*Wt]);
d2  = reshape(  d                         ,[N E* R*1*1*V*W ]); 
% d2  = reshape(  d - mean(d,1)             ,[N E* R*1*1*V*W ]); % removing the mean here


d2 = permute(d2,[2 1]);
J = d2*tp2; % [V E*R*K*F*1*W]
% [E*R*1*1*V*W Et*R*K*F*1*Wt]

J = reshape(J,[E*R*1*1*V*W Et R K F 1 Wt]);
J = permute(J,[2 3 4 5 6 7 1]);
J = reshape(J,[Et R K F 1 Wt E R 1 1 V W]);
J = permute(J,[13 7 2 3 4 11 6 1 5 8 9 10 12]); %[N ,E ,R ,K,F,V,W]
if any(size(J,9:15)~=1); dbstack; error('X'); end
% J = reshape(J,[V E R K F 1 W]); % [V E R K F 1 W]
% J = permute(J,[6 2 3 4 5 1 7]); % [N E R K F V W]