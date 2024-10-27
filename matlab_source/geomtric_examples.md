Euler Angles (Not Fixed) (eg. ZX'Z'')
%%% Direct 
phi = deg2rad(40)
theta = deg2rad(60)
psi = deg2rad(80)

R_zxz = euler_rotation('zxz', [phi, theta, psi])

%Equivalent to:
R = elem_rot_mat('z', phi)*elem_rot_mat('x', theta)*elem_rot_mat('z', psi)

%%% Inverse
[phi_inv, theta_inv, psi_inv] = euler_rotation_inverse('zxz', R_zxz, "pos")

fprintf("-----------------------------------------------------------------------")

RPY angles (Fixes XYZ)
%%% Direct
psi = deg2rad(80) % ROT X
theta = deg2rad(60) % ROT Y
phi = deg2rad(40) % ROT Z

R_rpy = rpy_rotation("xyz", [psi, theta, phi])

%Equivalent to:
R = elem_rot_mat('z', phi)*elem_rot_mat('y', theta)*elem_rot_mat('x', psi)

%%% Inverse
[phi_inv, theta_inv, psi_inv] = rpy_rotation_inverse("xyz", R_rpy, "pos")

angle_axis_orientation.mlx

r = [-1 1 1]'
r = r/norm(r)

theta = deg2rad(65)
R = angle_axis_rotation_direct(r, theta)

[theta_inv, r_inv] = get_theta_r(R, "pos")
theta_inv = rad2deg(theta_inv)

%singular case not managed
R = [-1 0 0; 
    0 -1/sqrt(2) -1/sqrt(2);
    0 -1/sqrt(2) 1/sqrt(2)]

[theta_inv, r_inv] = get_theta_r(R, "pos")

euler_rot_mat.mlx

%sequence of fixes axes ZYX
alpha_1 = -pi/2 %referred to Z
alpha_2 = -pi/4 %referred to Y
alpha_3 = pi/4  %referred to X

% R_X(alpha_3)*R_Y(alpha_2)*R_Z(alpha_1)
R_ZYX = euler_rotation('xyz', [alpha_3 alpha_2 alpha_1])

% how to get a unitary norm vector that will not be rotated by R
% i.e. R*r = r //rotation doesn't affect the vector
[V, D]=eig(R_ZYX) % V=eigenvector, D=eigenvalue
r = real(V(:,3)) %eigenvector associated to eigenvalue=1

euler_inverse.m

%sequence of fixes axes ZYX
alpha_1 = -pi/2 %referred to Z
alpha_2 = -pi/4 %referred to Y
alpha_3 = pi/4  %referred to X

% R_X(alpha_3)*R_Y(alpha_2)*R_Z(alpha_1)
R_ZYX = euler_rotation('xyz', [alpha_3 alpha_2 alpha_1])

[phi, theta, psi] = euler_rotation_inverse('xyz', R_ZYX, 'pos')