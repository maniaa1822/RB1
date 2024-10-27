# fuctions

function [T, A] = DHMatrix(arrays)
% T = DHMatrix(arrays) takes as inputs:
%   -arrays: a n-vector of vectors composed like this: [alpha a d theta]
% and outputs:
%   -T: the product of all the matrices corresponding to each vector of arrays
% Remember that:
% cos(q1 + q2) = cos(q1)*cos(q2) - sin(q1)*sin(q2)
% sin(q1 + q2) = cos(q1)*sin(q2) + cos(q2)*sin(q1)
% making use of the simplify function these are converted automatically

    T = eye(4);
    nums = size(arrays);
    
    A = cell(1,nums(1));
    
    for i = 1:nums(1)
        line = arrays(i, :);
        R = [cos(line(4)) -cos(line(1))*sin(line(4)) sin(line(1))*sin(line(4)) line(2)*cos(line(4));
             sin(line(4)) cos(line(1))*cos(line(4)) -sin(line(1))*cos(line(4)) line(2)*sin(line(4));
             0 sin(line(1)) cos(line(1)) line(3);
             0 0 0 1;];
        A{i} = R;
        T = T * R;   
    end

    if isa(T, 'sym')
        T = simplify(T);
    end
end

function f_r = get_f_r(T)
% f_r = get_f_r(T) takes as inputs:
%   -T: symbolic function of a homougeneous rototranslation matrix
%   -alpha_z: the fourth coordinate of the f_r mapping
%   -q: array of symbolic 
% and outputs:
%   -f_r: mapping from joint space to cartesian space

    if ~isa(T, 'sym')
        disp("We need a sym homogeneous matrix.");
        f_r = -1;
        return
    end
    
    alpha_z = sum(symvar(T(1:3, 1:3)));
    f_r = [T(1, 4); T(2, 4); T(3, 4); alpha_z];

end

% how to obtain f_r, (i.e. r)
%{
    From DH table:
        A0_2 = A0_1 * A1_2 //mapping from frame 2 to frame 0

    So, to obtain the origin of frame 2 in frame 0:
        O2_2 = [0 0 0] % 3D space
            -> O2_2_H = [0 0 0 1] //homog. transf
        
        r = O0_2 = [rx ry rz 1]' = A0_2*O2_2_H 
%}

## use example of get_f_r() and DHMatrix()

syms q1 real
syms q2 real
syms q3 real 

syms l1 real

alpha = [0 pi/2 0];
a=[l1 0 0];
d=[0 0 q3];
theta=[q1 q2+pi/2 0];

table=[alpha',a',d',theta']

[T, A] = DHMatrix(table);
A0_1=A{1}
A1_2=A{2}
A2_3=A{3}

T
% how to obtain f_r, (i.e. r)
%{
    From DH table:
        A0_2 = A0_1 * A1_2 //mapping from frame 2 to frame 0
    So, to obtain the origin of frame 2 in frame 0:
        O2_2 = [0 0 0] % 3D space
            -> O2_2_H = [0 0 0 1] //homog. transf
        
        r = O0_2 = [rx ry rz 1]' = A0_2*O2_2_H 
%}
f_r_3D = get_f_r(T); %(X Y Z PHI)

f_r = [f_r_3D(1); f_r_3D(2); f_r_3D(4)]

j = jacobian(f_r, [q1 q2 q3])
j0_3 = j;
j1_3 = (get_rotation_mat(A0_1))' * j0_3;
j1_3 = simplify(j1_3)

simplify(det(j1_3))

%{
minors = get_minors(j1_3, 2)
det_1 = simplify(det(simplify(minors{1})))
det_2 = simplify(det(simplify(minors{2})))
det_3 = simplify(det(simplify(minors{3})))
det_4 = simplify(det(simplify(minors{4})))
det_5 = simplify(det(simplify(minors{5})))
det_6 = simplify(det(simplify(minors{6})))
det_7 = simplify(det(simplify(minors{7})))
det_8 = simplify(det(simplify(minors{8})))
det_9 = simplify(det(simplify(minors{9})))
%}

f = [0 1.5 -4.5]'
j_0 = subs(j, {q1, q2, q3, l1}, {pi/2, 0, 3, 0.5})

tau_0 = j_0*f
tau_0 = -tau_0


f = [0 1.5 -4.5]'
j_0 = subs(j, {q1, q2, q3, l1}, {0, pi/2,0, 0.5})

tau_s = j_0*f
tau_s = -tau_s