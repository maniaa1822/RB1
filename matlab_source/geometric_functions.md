function R = rpy_rotation(sequence, angles)
% R = rpy_rotation(sequence, angles) takes as inputs:
%   -sequence: a string which specifies which elementary rotation matrices
%              must be multiplied together to obtain the desired rotation
%   -angles: the radiants of the three rotation angles
% and outputs:
%   -R: The desired rotation
% RPY rotations work about fixed-axes

    sequence = char(sequence);
    R = euler_rotation(flip(sequence), flip(angles));
end

function [phi, theta, psi] = rpy_rotation_inverse(sequence, R, pos_or_neg_solution)
% [phi, theta, psi] = rpy_rotation_inverse(sequence, R) takes as inputs:
%   -sequence: a string which specifies how the RPY-rotation has been computed, e.g. "xyx"
%   -R: the rotation to be decomposed, should be a 3x3 matrix
% and outputs:
%   -phi: the radiants of the first rotation
%   -theta: the radiants of the second rotation
%   -psi: the radiants of the third rotation
% Upon execution the function will ask for an input which will determine
% which solution must be displayed (there is always a pair of solutions)
%
% RPY rotations work about fixed-axes

    sequence = char(sequence);
    %[phi, theta, psi] = euler_rotation_inverse(flip(sequence), R, pos_or_neg_solution);
    [psi, theta, phi] = euler_rotation_inverse(flip(sequence), R, pos_or_neg_solution);

function S_s = M_to_S_skew_sym(M, q, q_dot, print_info)
    cell_c_k = christoffel_symbols(M, q);

    if print_info
        fprintf("----------------------------------------------------\n");
        fprintf("Matrices of Christoffel symbols (c_i):\n");
        print_cells(cell_c_k)
        fprintf("\n");
    end
    
    S_s = sym([]);
    
    n = length(q); 
    for i=1:n
        s_i = q_dot'*cell_c_k{i};
        S_s = [S_s; s_i];
    end
    
    if print_info
        fprintf("Row i-th of S matrix is given by s_i = q_dot' * c_i\n");
        fprintf("Result S:\n");
        display(S_s);
        fprintf("----------------------------------------------------\n");
    end
end

function C = M_to_C(M, q, q_dot)
    C = sym([]);
    n = length(q);
    
    for i=1:n
        M_k = M(:, i);
        c_k = 0.5*(jacobian(M_k, q)+(jacobian(M_k, q)') - diff(M, q(i)));
        c_k = simplify(c_k);
        C = [C; simplify((q_dot')*c_k*q_dot)];
    end
end

function res = isRotationMatrix(R)
    RR = eval(simplify(R'*R));
    det_R = eval(simplify(det(R)));
    if isequal(RR, eye(3)) && isequal(det_R, 1)
        res = true;
    else
        res = false;
    end
end

    function R = euler_rotation(sequence, angles)
% R = euler_rotation(sequence, angles) takes as inputs:
%   -sequence: a string which specifies the axes along which rotation
%              occurs
%   -angles: The radiants (or symbolics) of the three rotations
% and outputs:
%   -R: The desired rotation
% Euler rotations work about moving-axes
    
    if strlength(sequence) ~= 3
        disp("Sequence not valid, must be of lenght three.")
        return;
    end
    
    sequence = lower(char(sequence));
    if (sequence(2) == sequence(1) || sequence(2) == sequence(3))
        disp("Two consecutive rotation along the same axis are not valid.")
        return
    end
    
    R = elem_rot_mat(sequence(1), angles(1)) * elem_rot_mat(sequence(2), angles(2)) * elem_rot_mat(sequence(3), angles(3));

end

function [phi, theta, psi] = euler_rotation_inverse(sequence, R, pos_or_neg_solution)
% [phi, theta, psi] = euler_rotation_inverse(sequence, R) takes as inputs:
%   -sequence: a string which specifies how the euler-rotation has been computed, e.g. "xyx"
%   -R: the rotation to be decomposed, should be a 3x3 matrix
% and outputs:
%   -phi: the radiants of the first rotation
%   -theta: the radiants of the second rotation
%   -psi: the radiants of the third rotation
% Upon execution the function will ask for an input which will determine
% which solution must be displayed (there is always a pair of solutions)
%
% Euler rotations work about moving-axes
   
%    if isRotationMatrix(R) == false
%        fprintf("ERROR: R is not a rotation matrix!")
%        return
%    end
    
    if isa(R, 'sym')
        disp("R must by non-symbolic matrix")
        phi = -1; theta = -1; psi = -1;
        return
    end
    
    if strlength(sequence) ~= 3
        disp("Sequence not valid, must be of length three.")
        return;
    end
    
    if (sequence(2) == sequence(1) || sequence(2) == sequence(3))
        disp("Two consecutive rotation along the same axis are not valid.")
        return
    end
    
    risp = pos_or_neg_solution; %input("There is always a pair of solutions to the inverse euler problem. Do you want to use the positive or negative sin(theta)? (pos, neg)\n", "s");
    cond = risp == "pos";
    
    switch lower(sequence)
        case "xyx"
            theta = atan2(sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*cond + atan2(-sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(1, 2)/sin(theta), R(1, 3)/sin(theta));
            phi = atan2(R(2, 1)/sin(theta), -R(3, 1)/sin(theta));
            
        case "xyz"
            theta = atan2(R(1, 3), sqrt(R(1, 1)^2 + R(1, 2)^2))*cond + atan2(R(1, 3), -sqrt(R(1, 1)^2 + R(1, 2)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(-R(1, 2)/cos(theta), R(1, 1)/cos(theta));
            phi = atan2(-R(2, 3)/cos(theta), R(3, 3)/cos(theta));
            
        case "xzx"
            theta = atan2(sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*cond + atan2(-sqrt(R(1, 2)^2 + R(1, 3)^2), R(1, 1))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(1, 3)/sin(theta), -R(1, 2)/sin(theta));
            phi = atan2(R(3, 1)/sin(theta), R(2, 1)/sin(theta));
            
        case "xzy"
            theta = atan2(-R(1, 2), sqrt(R(1, 1)^2 + R(1, 3)^2))*cond + atan2(-R(1, 1), -sqrt(R(1, 3)^2 + R(2, 1)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(1, 3)/cos(theta), R(1, 1)/cos(theta));
            phi = atan2(R(3, 2)/cos(theta), R(2, 2)/cos(theta));
            
        case "yxy"
            theta = atan2(sqrt(R(2, 3)^2 + R(2, 1)^2), R(2, 2))*cond + atan2(-sqrt(R(2, 3)^2 + R(2, 1)^2), R(2, 2))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(2, 1)/sin(theta), -R(2, 3)/sin(theta));
            phi = atan2(R(1, 2)/sin(theta), R(3, 2)/sin(theta));
            
        case "yxz"
            theta = atan2(-R(2, 3), sqrt(R(2, 2)^2 + R(2, 1)^2))*cond + atan2(-R(2, 3), -sqrt(R(2, 2)^2 + R(2, 1)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(2, 1)/cos(theta), R(2, 2)/cos(theta));
            phi = atan2(R(1, 3)/cos(theta), R(3, 3)/cos(theta));
            
        case "yzx"
            theta = atan2(R(2, 1), sqrt(R(2, 2)^2 + R(2, 3)^2))*cond + atan2(R(2, 1), -sqrt(R(2, 2)^2 + R(2, 3)^2))*(1-cond);

            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(-R(2, 3)/cos(theta), R(2, 2)/cos(theta));
            phi = atan2(-R(3, 1)/cos(theta), R(1, 1)/cos(theta));
            
        case "yzy"
            theta = atan2(sqrt(R(2, 1)^2 + R(2,3)^2), R(2, 2))*cond + atan2(-sqrt(R(2, 1)^2 + R(2,3)^2), R(2, 2))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(2, 3)/sin(theta), R(2, 1)/sin(theta));
            phi = atan2(R(3, 2)/sin(theta), -R(1, 2)/sin(theta));
            
        case "zxy"
            theta = atan2(R(3, 2), sqrt(R(3, 1)^2 + R(3,3)^2))*cond + atan2(R(3, 2), -sqrt(R(3, 1)^2 + R(3,3)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(-R(3, 1)/cos(theta), R(3, 3)/cos(theta));
            phi = atan2(-R(1, 2)/cos(theta), R(2, 1)/cos(theta));
            
        case "zxz"
            theta = atan2(sqrt(R(1, 3)^2 + R(2, 3)^2), R(3, 3))*cond + atan2(-sqrt(R(1, 3)^2 + R(2, 3)^2), R(3, 3))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(3, 1)/sin(theta), R(3, 2)/sin(theta));
            phi = atan2(R(1, 3)/sin(theta), -R(2, 3)/sin(theta));

        case "zyx"
            theta = atan2(-R(3, 1), sqrt(R(3, 2)^2+R(3, 3)^2))*cond + atan2(-R(3, 1), -sqrt(R(3, 2)^2+R(3, 3)^2))*(1-cond);
            if (abs(cos(theta)) <= 1e-6)
                disp('Singular case: cos(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(3, 2)/cos(theta), R(3, 3)/cos(theta));
            phi = atan2(R(2, 1)/cos(theta), R(1, 1)/cos(theta));
            
        case "zyz"
            theta = atan2(sqrt(R(3, 1)^2 + R(3, 2)^2), R(3, 3))*cond + atan2(-sqrt(R(3, 1)^2 + R(3, 2)^2), R(3, 3))*(1-cond);
            if (abs(sin(theta)) <= 1e-6)
                disp('Singular case: sin(theta) == 0 or very close to 0.')
                return
            end
            
            psi = atan2(R(3, 2)/sin(theta), -R(3, 1)/sin(theta));
            phi = atan2(R(2, 3)/sin(theta), R(1, 3)/sin(theta));
            
        otherwise
            disp("Invalid sequence")
    end
    
end

function R = elem_rot_mat(axis, s)
% R = elem_rot_mat(axis, angle) takes as inputs:
%   -axis: The axis we should perform the rotation about, ("x", "y", "z")
%   -angle: The radiants (or a symbolic) we should perform the rotation of
% and outputs:
%   -R: The desired elementary rotation

    switch axis
        case {"x", "X"}
            R = [1       0        0;
                 0    cos(s)   -sin(s);
                 0    sin(s)    cos(s)];
        case {"y", "Y"}
            R = [cos(s)     0   sin(s);
                   0	    1     0;
                -sin(s)     0   cos(s)];
        case {"z", "Z"}
            R = [cos(s)  -sin(s)    0;
                 sin(s)   cos(s)    0;
                   0        0       1];
        otherwise
            disp("First Parameter should be either 'x', 'y', 'z' or any of those capitalized")
    end
    
end

function alpha = elem_rot_mat_inverse(axis, rotation_matrix)
    switch axis
        case {"x", "X"}
%             R = [1       0        0;
%                  0    cos(s)   -sin(s);
%                  0    sin(s)    cos(s)];
             
            sin_alpha = rotation_matrix(3, 2);
            cos_alpha = rotation_matrix(3, 3);
            alpha = atan2(sin_alpha, cos_alpha);
            
        case {"y", "Y"}            
%             R = [cos(s)     0   sin(s);
%                     0	    1     0;
%                 -sin(s)     0   cos(s)];

            sin_alpha = rotation_matrix(1, 3);
            cos_alpha = rotation_matrix(1, 1);
            alpha = atan2(sin_alpha, cos_alpha);
        
        case {"z", "Z"}
%             R = [cos(s)  -sin(s)    0;
%                  sin(s)   cos(s)    0;
%                    0        0       1];

            sin_alpha = rotation_matrix(2, 1);
            cos_alpha = rotation_matrix(1, 1);
            alpha = atan2(sin_alpha, cos_alpha);
        otherwise
            disp("First Parameter should be either 'x', 'y', 'z' or any of those capitalized")
    end

end

function R = angle_axis_rotation_direct(r, theta)
% angle-axis rotation representation, DIRECT PROBLEM 
    I = eye(3);
    R = r*r' + (I-r*r')*cos(theta)+S(r)*sin(theta);
end

function T = affine_T(rotation_matrix, translation)
    T = [rotation_matrix, translation; 0 0 0 1];
end

function translation = affine_get_translation(T)
    translation = [T(1, 4); T(2, 4); T(3, 4)];
end

function R = affine_get_R(T)
    R = T(1:3, 1:3);
end
