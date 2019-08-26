function p_transformed = Transform_3D_matrix(Tr,p)


 
[N_puntos,aux] = size(p); 


p_transformed = zeros(size(p));



p_Homogeneous = [[p],ones(N_puntos,1)];
p_tr_Homogeneous = p_Homogeneous;
p_tr_Homogeneous(:,1:3)=0.0;

p_tr_Homogeneous = Tr*p_Homogeneous';

p_transformed(:,:)=(p_tr_Homogeneous(1:3,:))';


end