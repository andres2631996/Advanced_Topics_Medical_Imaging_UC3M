function image_new = transform_image_3D (Tr,image, interp_type)

[M,N,L] = size(image);


% Inverse Matrix necessary for the interpolation
Tr_inv = inv(Tr);

[X,Y,Z] = meshgrid(1:N,1:M,1:L);


points_transf_inv = Transform_3D_matrix(Tr_inv,[X(:),Y(:),Z(:)]);

X_transf_inv = reshape(points_transf_inv (:,1),size(X));
Y_transf_inv = reshape(points_transf_inv (:,2),size(Y));
Z_transf_inv = reshape(points_transf_inv (:,3),size(Z));


% Interpolation is performed
image_new = interp3(image,X_transf_inv,Y_transf_inv,Z_transf_inv,interp_type);
image_new(isnan(image_new)) = 0;


end