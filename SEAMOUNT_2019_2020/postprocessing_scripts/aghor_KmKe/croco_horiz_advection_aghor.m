function Tadv=croco_horiz_advection_aghor(masku,maskv,maskr,Hz,dn_uXu,dm_vXv,temp)
% Aghor: only horizontal advection, no vertical term here

%
% Compute Uflx: 
%
FlxU=0.5*(Hz(:,:,2:end)+Hz(:,:,1:end-1)).*dn_uXu;
%
% Compute Vflx
%
FlxV=0.5*(Hz(:,2:end,:)+Hz(:,1:end-1,:)).*dm_vXv;
%
% Compute Flux term and "diffusion" in xi direction
%
[Xadv,trun]=croco_horiz_trcflux(masku,maskr,temp,FlxU,1);
%
% Compute Flux term and "diffusion" in eta direction
%
[Yadv,trun]=croco_horiz_trcflux(maskv,maskr,temp,FlxV,2);
%
Tadv=-Xadv-Yadv;
%
return
