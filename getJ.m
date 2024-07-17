function J = getJ(data,tp,tvec,f)
% data [time x vox   x trial]
% tp   [time x trial x taper]
% tvec [time x trial x taper]
% f    [1    x 1     x freq ]
[N,C,E,R]  = size(data);
[Nk,Ek,K]  = size(tp);
[Nv,Ev,Kv]  = size(tvec);
F = length(f);
if Kv==1 && Kv~=K
    tvec = repmat(tvec,[1 1 K]);
end
clear Nk Ek Nv Ev Kv

tvec2 = exp(  f  .*  reshape(tvec*2*pi*1i,[N*E E*K])  ); %[N*E x E*K   x F]
tapers2 = reshape(  reshape(tp,[N*E E*K])  .*  tvec2  ,[N*E E*K*F]);  %[N*E x E*K*F    ]
data2 = reshape(permute(data,[2 1 3]),[C N*E]);
J = permute(  reshape(  data2 * tapers2  ,[C E K F])  ,[4 3 1 2]); % [freq x taper x vox x trial]
