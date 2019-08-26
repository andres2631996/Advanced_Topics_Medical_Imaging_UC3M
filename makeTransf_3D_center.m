function Tr = makeTransf_3D_center(rx,ry,rz,tx,ty,tz,sx,sy,sz,center)


Rz = [cosd(rz) sind(rz) 0 0; -sind(rz) cosd(rz) 0 0; 0 0 1 0; 0 0 0 1];
Rx = [ 1 0 0 0; 0 cosd(rx) sind(rx) 0; 0 -sind(rx) cosd(rx) 0; 0 0 0 1];
Ry = [cosd(ry) 0 sind(ry) 0; 0 1 0 0; -sind(ry) 0 cosd(ry) 0; 0 0 0 1];
%R = [cosd(angle) -sind(angle) 0; sind(angle) cosd(angle) 0; 0 0 1];
T = [1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];
S = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];

T_center = [1 0 0 -center(2); 0 1 0 -center(1); 0 0 1 -center(3); 0 0 0 1];
%T_center = [1 0 -center(1); 0 1 -center(2); 0 0 1];
T_center_inv = inv(T_center);
%T_center_inv = [1 0 center(1); 0 1 center(2); 0 0 1];

Tr = T_center_inv*S*Rx*Ry*Rz*T*T_center;


end