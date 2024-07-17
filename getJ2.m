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


tp2 = reshape(  tp .* exp(-f.*t*2*pi*1i)  ,[N Et*R*K*F*1*Wt]);
d2  = reshape(  d - mean(d,1)             ,[N  E*R*1*1*V*W ]); % removing the mean here

d2 = permute(d2,[2 1]);
J = d2*tp2; % [V E*R*K*F*1*W]

J = reshape(J,[V E R K F 1 W]); % [V E R K F 1 W]
J = permute(J,[6 2 3 4 5 1 7]); % [N E R K F V W]