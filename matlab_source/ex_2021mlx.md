
```matlab
%rotation with Eulero angles (moving axes)
%so R_YXZ = R_Y*R_X*R_Z
syms alpha real
syms beta real
syms gamma real

R_yxz = euler_rotation('yxz', [alpha beta gamma])
```
R_yxz = 
 $\displaystyle \left(\begin{array}{ccc} \cos \left(\alpha \right)\,\cos \left(\gamma \right)+\sin \left(\alpha \right)\,\sin \left(\beta \right)\,\sin \left(\gamma \right) & \cos \left(\gamma \right)\,\sin \left(\alpha \right)\,\sin \left(\beta \right)-\cos \left(\alpha \right)\,\sin \left(\gamma \right) & \cos \left(\beta \right)\,\sin \left(\alpha \right)\newline \cos \left(\beta \right)\,\sin \left(\gamma \right) & \cos \left(\beta \right)\,\cos \left(\gamma \right) & -\sin \left(\beta \right)\newline \cos \left(\alpha \right)\,\sin \left(\beta \right)\,\sin \left(\gamma \right)-\cos \left(\gamma \right)\,\sin \left(\alpha \right) & \sin \left(\alpha \right)\,\sin \left(\gamma \right)+\cos \left(\alpha \right)\,\cos \left(\gamma \right)\,\sin \left(\beta \right) & \cos \left(\alpha \right)\,\cos \left(\beta \right) \end{array}\right)$
 

```matlab
b = atan2(-sqrt(3)/2, 0.5)
```

```matlabTextOutput
b = -1.0472
```

```matlab
a = atan2(0, (-0.5)/cos(b))
```

```matlabTextOutput
a = 3.1416
```

```matlab
c = atan2(0.5/cos(b), 0)
```

```matlabTextOutput
c = 1.5708
```

```matlab

a = rad2deg(a)
```

```matlabTextOutput
a = 180
```

```matlab
b=rad2deg(b)
```

```matlabTextOutput
b = -60.0000
```

```matlab
c=rad2deg(c)
```

```matlabTextOutput
c = 90
```

```matlab
R = subs(R_yxz, {alpha, beta, gamma}, {deg2rad(180), deg2rad(-60), deg2rad(90)})
```
R = 
 $\displaystyle \left(\begin{array}{ccc} 0 & 1 & 0\newline \frac{1}{2} & 0 & \frac{\sqrt{3}}{2}\newline \frac{\sqrt{3}}{2} & 0 & -\frac{1}{2} \end{array}\right)$
 

```matlab
WR_B = R
```
WR_B = 
 $\displaystyle \left(\begin{array}{ccc} 0 & 1 & 0\newline \frac{1}{2} & 0 & \frac{\sqrt{3}}{2}\newline \frac{\sqrt{3}}{2} & 0 & -\frac{1}{2} \end{array}\right)$
 

```matlab
BV_B = [1 -1 0]'
```

```matlabTextOutput
BV_B = 3x1
     1
    -1
     0

```

```matlab
WR_B_dot = S(WR_B*BV_B)*WR_B
```
WR_B_dot = 
 $\displaystyle \left(\begin{array}{ccc} 0 & 0 & -1\newline \frac{\sqrt{3}}{2} & \frac{\sqrt{3}}{2} & -\frac{1}{2}\newline -\frac{1}{2} & -\frac{1}{2} & -\frac{\sqrt{3}}{2} \end{array}\right)$
 

```matlab
syms q1 real
syms q2 real
syms q3 real 

syms l0 real
syms l1 real
syms l2 real
syms l3 real

alpha = [-pi/2 -pi/2 0];
a=[-l1 -l2 -l3];
d=[l0 0 0];
theta=[q1 q2 q3];

table=[alpha',a',d',theta']
```
table = 
 $\displaystyle \left(\begin{array}{cccc} -\frac{\pi }{2} & -l_1  & l_0  & q_1 \newline -\frac{\pi }{2} & -l_2  & 0 & q_2 \newline 0 & -l_3  & 0 & q_3  \end{array}\right)$
 

```matlab

[T, A] = DHMatrix(table);
A0_1=A{1}
```
A0_1 = 
 $\displaystyle \left(\begin{array}{cccc} \cos \left(q_1 \right) & 0 & -\sin \left(q_1 \right) & -l_1 \,\cos \left(q_1 \right)\newline \sin \left(q_1 \right) & 0 & \cos \left(q_1 \right) & -l_1 \,\sin \left(q_1 \right)\newline 0 & -1 & 0 & l_0 \newline 0 & 0 & 0 & 1 \end{array}\right)$
 

```matlab
A1_2=A{2}
```
A1_2 = 
 $\displaystyle \left(\begin{array}{cccc} \cos \left(q_2 \right) & 0 & -\sin \left(q_2 \right) & -l_2 \,\cos \left(q_2 \right)\newline \sin \left(q_2 \right) & 0 & \cos \left(q_2 \right) & -l_2 \,\sin \left(q_2 \right)\newline 0 & -1 & 0 & 0\newline 0 & 0 & 0 & 1 \end{array}\right)$
 

```matlab
A2_2=A{3}
```
A2_2 = 
 $\displaystyle \left(\begin{array}{cccc} \cos \left(q_3 \right) & -\sin \left(q_3 \right) & 0 & -l_3 \,\cos \left(q_3 \right)\newline \sin \left(q_3 \right) & \cos \left(q_3 \right) & 0 & -l_3 \,\sin \left(q_3 \right)\newline 0 & 0 & 1 & 0\newline 0 & 0 & 0 & 1 \end{array}\right)$
 

```matlab
T
```
T = 
 $\displaystyle \left(\begin{array}{cccc} \sin \left(q_1 \right)\,\sin \left(q_3 \right)+\cos \left(q_1 \right)\,\cos \left(q_2 \right)\,\cos \left(q_3 \right) & \cos \left(q_3 \right)\,\sin \left(q_1 \right)-\cos \left(q_1 \right)\,\cos \left(q_2 \right)\,\sin \left(q_3 \right) & -\cos \left(q_1 \right)\,\sin \left(q_2 \right) & -l_1 \,\cos \left(q_1 \right)-l_2 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)-l_3 \,\sin \left(q_1 \right)\,\sin \left(q_3 \right)-l_3 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\newline \cos \left(q_2 \right)\,\cos \left(q_3 \right)\,\sin \left(q_1 \right)-\cos \left(q_1 \right)\,\sin \left(q_3 \right) & -\cos \left(q_1 \right)\,\cos \left(q_3 \right)-\cos \left(q_2 \right)\,\sin \left(q_1 \right)\,\sin \left(q_3 \right) & -\sin \left(q_1 \right)\,\sin \left(q_2 \right) & l_3 \,\cos \left(q_1 \right)\,\sin \left(q_3 \right)-l_2 \,\cos \left(q_2 \right)\,\sin \left(q_1 \right)-l_1 \,\sin \left(q_1 \right)-l_3 \,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\,\sin \left(q_1 \right)\newline -\cos \left(q_3 \right)\,\sin \left(q_2 \right) & \sin \left(q_2 \right)\,\sin \left(q_3 \right) & -\cos \left(q_2 \right) & l_0 +l_2 \,\sin \left(q_2 \right)+l_3 \,\cos \left(q_3 \right)\,\sin \left(q_2 \right)\newline 0 & 0 & 0 & 1 \end{array}\right)$
 

```matlab
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

f_r = [f_r_3D(1); f_r_3D(2); f_r_3D(3)]
```
f_r = 
 $\displaystyle \left(\begin{array}{c} -l_1 \,\cos \left(q_1 \right)-l_2 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)-l_3 \,\sin \left(q_1 \right)\,\sin \left(q_3 \right)-l_3 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\newline l_3 \,\cos \left(q_1 \right)\,\sin \left(q_3 \right)-l_2 \,\cos \left(q_2 \right)\,\sin \left(q_1 \right)-l_1 \,\sin \left(q_1 \right)-l_3 \,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\,\sin \left(q_1 \right)\newline l_0 +l_2 \,\sin \left(q_2 \right)+l_3 \,\cos \left(q_3 \right)\,\sin \left(q_2 \right) \end{array}\right)$
 

```matlab
op = f_r 
```
op = 
 $\displaystyle \left(\begin{array}{c} -l_1 \,\cos \left(q_1 \right)-l_2 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)-l_3 \,\sin \left(q_1 \right)\,\sin \left(q_3 \right)-l_3 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\newline l_3 \,\cos \left(q_1 \right)\,\sin \left(q_3 \right)-l_2 \,\cos \left(q_2 \right)\,\sin \left(q_1 \right)-l_1 \,\sin \left(q_1 \right)-l_3 \,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\,\sin \left(q_1 \right)\newline l_0 +l_2 \,\sin \left(q_2 \right)+l_3 \,\cos \left(q_3 \right)\,\sin \left(q_2 \right) \end{array}\right)$
 

```matlab
wRo = [-1 0 0;
        0 0 1;
        0 1 0]
```

```matlabTextOutput
wRo = 3x3
    -1     0     0
     0     0     1
     0     1     0

```

```matlab
wp = wRo * op 
```
wp = 
 $\displaystyle \left(\begin{array}{c} l_1 \,\cos \left(q_1 \right)+l_2 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)+l_3 \,\sin \left(q_1 \right)\,\sin \left(q_3 \right)+l_3 \,\cos \left(q_1 \right)\,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\newline l_0 +l_2 \,\sin \left(q_2 \right)+l_3 \,\cos \left(q_3 \right)\,\sin \left(q_2 \right)\newline l_3 \,\cos \left(q_1 \right)\,\sin \left(q_3 \right)-l_2 \,\cos \left(q_2 \right)\,\sin \left(q_1 \right)-l_1 \,\sin \left(q_1 \right)-l_3 \,\cos \left(q_2 \right)\,\cos \left(q_3 \right)\,\sin \left(q_1 \right) \end{array}\right)$
 

```matlab
syms q1 real
syms q2 real
syms q3 real 
syms q4 real 

syms l1 real
syms l2 real

alpha = [0 -pi/2 pi/2 0];
a=[l1 0 0 l2];
d=[0 0 q3 0];
theta=[q1 q2-pi/2 0 q4+pi/2];

table=[alpha',a',d',theta']
```
table = 
 $\displaystyle \left(\begin{array}{cccc} 0 & l_1  & 0 & q_1 \newline -\frac{\pi }{2} & 0 & 0 & q_2 -\frac{\pi }{2}\newline \frac{\pi }{2} & 0 & q_3  & 0\newline 0 & l_2  & 0 & q_4 +\frac{\pi }{2} \end{array}\right)$
 

```matlab

[T, A] = DHMatrix(table);
A0_1=A{1}
```
A0_1 = 
 $\displaystyle \left(\begin{array}{cccc} \cos \left(q_1 \right) & -\sin \left(q_1 \right) & 0 & l_1 \,\cos \left(q_1 \right)\newline \sin \left(q_1 \right) & \cos \left(q_1 \right) & 0 & l_1 \,\sin \left(q_1 \right)\newline 0 & 0 & 1 & 0\newline 0 & 0 & 0 & 1 \end{array}\right)$
 

```matlab
A1_2=A{2}
```
A1_2 = 
 $\displaystyle \left(\begin{array}{cccc} \cos \left(q_2 -\frac{\pi }{2}\right) & 0 & -\sin \left(q_2 -\frac{\pi }{2}\right) & 0\newline \sin \left(q_2 -\frac{\pi }{2}\right) & 0 & \cos \left(q_2 -\frac{\pi }{2}\right) & 0\newline 0 & -1 & 0 & 0\newline 0 & 0 & 0 & 1 \end{array}\right)$
 

```matlab
A2_2=A{3}
```
A2_2 = 
 $\displaystyle \left(\begin{array}{cccc} 1 & 0 & 0 & 0\newline 0 & 0 & -1 & 0\newline 0 & 1 & 0 & q_3 \newline 0 & 0 & 0 & 1 \end{array}\right)$
 

```matlab
T
```
T = 
 $\displaystyle \begin{array}{l} \left(\begin{array}{cccc} \sigma_2  & -\sigma_1  & 0 & q_3 \,\cos \left(q_1 +q_2 \right)+l_1 \,\cos \left(q_1 \right)+l_2 \,\sigma_2 \newline \sigma_1  & \sigma_2  & 0 & q_3 \,\sin \left(q_1 +q_2 \right)+l_1 \,\sin \left(q_1 \right)+l_2 \,\sigma_1 \newline 0 & 0 & 1 & 0\newline 0 & 0 & 0 & 1 \end{array}\right)\\\mathrm{}\\\textrm{where}\\\mathrm{}\\\;\;\sigma_1 =\sin \left(q_1 +q_2 +q_4 \right)\\\mathrm{}\\\;\;\sigma_2 =\cos \left(q_1 +q_2 +q_4 \right)\end{array}$
 

```matlab
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
```
f_r = 
 $\displaystyle \left(\begin{array}{c} q_3 \,\cos \left(q_1 +q_2 \right)+l_1 \,\cos \left(q_1 \right)+l_2 \,\cos \left(q_1 +q_2 +q_4 \right)\newline q_3 \,\sin \left(q_1 +q_2 \right)+l_1 \,\sin \left(q_1 \right)+l_2 \,\sin \left(q_1 +q_2 +q_4 \right)\newline q_1 +q_2 +q_4  \end{array}\right)$
 

```matlab
j = simplify(jacobian(f_r, [q1 q2 q3 q4]))
```
j = 
 $\displaystyle \begin{array}{l} \left(\begin{array}{cccc} -q_3 \,\sin \left(q_1 +q_2 \right)-l_1 \,\sin \left(q_1 \right)-\sigma_1  & -q_3 \,\sin \left(q_1 +q_2 \right)-\sigma_1  & \cos \left(q_1 +q_2 \right) & -\sigma_1 \newline q_3 \,\cos \left(q_1 +q_2 \right)+l_1 \,\cos \left(q_1 \right)+\sigma_2  & q_3 \,\cos \left(q_1 +q_2 \right)+\sigma_2  & \sin \left(q_1 +q_2 \right) & \sigma_2 \newline 1 & 1 & 0 & 1 \end{array}\right)\\\mathrm{}\\\textrm{where}\\\mathrm{}\\\;\;\sigma_1 =l_2 \,\sin \left(q_1 +q_2 +q_4 \right)\\\mathrm{}\\\;\;\sigma_2 =l_2 \,\cos \left(q_1 +q_2 +q_4 \right)\end{array}$
 

```matlab
j_ref_1 = get_rotation_mat(A0_1)' * j
```
j_ref_1 = 
 $\displaystyle \begin{array}{l} \left(\begin{array}{cccc} \sin \left(q_1 \right)\,\sigma_2 -\cos \left(q_1 \right)\,\sigma_1  & \sin \left(q_1 \right)\,{\left(\sigma_6 +l_2 \,\sigma_5 \right)}-\cos \left(q_1 \right)\,{\left(\sigma_4 +l_2 \,\sigma_3 \right)} & \cos \left(q_1 +q_2 \right)\,\cos \left(q_1 \right)+\sin \left(q_1 +q_2 \right)\,\sin \left(q_1 \right) & l_2 \,\sigma_5 \,\sin \left(q_1 \right)-l_2 \,\sigma_3 \,\cos \left(q_1 \right)\newline \cos \left(q_1 \right)\,\sigma_2 +\sin \left(q_1 \right)\,\sigma_1  & \sin \left(q_1 \right)\,{\left(\sigma_4 +l_2 \,\sigma_3 \right)}+\cos \left(q_1 \right)\,{\left(\sigma_6 +l_2 \,\sigma_5 \right)} & \sin \left(q_1 +q_2 \right)\,\cos \left(q_1 \right)-\cos \left(q_1 +q_2 \right)\,\sin \left(q_1 \right) & l_2 \,\sigma_5 \,\cos \left(q_1 \right)+l_2 \,\sigma_3 \,\sin \left(q_1 \right)\newline 1 & 1 & 0 & 1 \end{array}\right)\\\mathrm{}\\\textrm{where}\\\mathrm{}\\\;\;\sigma_1 =\sigma_4 +l_1 \,\sin \left(q_1 \right)+l_2 \,\sigma_3 \\\mathrm{}\\\;\;\sigma_2 =\sigma_6 +l_1 \,\cos \left(q_1 \right)+l_2 \,\sigma_5 \\\mathrm{}\\\;\;\sigma_3 =\sin \left(q_1 +q_2 +q_4 \right)\\\mathrm{}\\\;\;\sigma_4 =q_3 \,\sin \left(q_1 +q_2 \right)\\\mathrm{}\\\;\;\sigma_5 =\cos \left(q_1 +q_2 +q_4 \right)\\\mathrm{}\\\;\;\sigma_6 =q_3 \,\cos \left(q_1 +q_2 \right)\end{array}$
 

```matlab
minors = get_minors(j_ref_1, 3)
```
| |1|2|3|4|
|:--:|:--:|:--:|:--:|:--:|
|1|3x3 sym|3x3 sym|3x3 sym|3x3 sym|

```matlab
det_1 = det(simplify(minors{1}))
```
det_1 = 
 $\displaystyle l_1 \,\cos \left(q_2 \right)$
 

```matlab
det_2 = det(simplify(minors{2}))
```
det_2 = 
 $\displaystyle l_1 \,q_3 \,\sin \left(q_2 \right)$
 

```matlab
det_3 = det(simplify(minors{3}))
```
det_3 = 
 $\displaystyle -q_3 \,{\cos \left(q_2 \right)}^2 -l_1 \,\cos \left(q_2 \right)-q_3 \,{\sin \left(q_2 \right)}^2 $
 

```matlab
det_4 = det(simplify(minors{4}))
```
det_4 = 
 $\displaystyle -q_3 \,{\cos \left(q_2 \right)}^2 -q_3 \,{\sin \left(q_2 \right)}^2 $
 

```matlab
rank_j=rank(j_ref_1)
```

```matlabTextOutput
rank_j = 3
```

```matlab
% Find the singular configurations (where rank(J) < min(rows, columns))
```

```matlabTextOutput
singular_configurations = struct with fields:
    q1: [0x1 sym]
    q2: [0x1 sym]
    q3: [0x1 sym]
    q4: [0x1 sym]

```

```matlab
singular_1 = solve(det_1 == 0, [q1 q2 q3 q4])
```

```matlabTextOutput
singular_1 = struct with fields:
    q1: 0
    q2: pi/2
    q3: 0
    q4: 0

```

```matlab
singular_2 = solve(det_2 == 0, [q1 q2 q3 q4])
```

```matlabTextOutput
singular_2 = struct with fields:
    q1: 0
    q2: 0
    q3: 0
    q4: 0

```

```matlab
singular_3 = solve(det_3 == 0, [q1 q2 q3 q4])
```

```matlabTextOutput
singular_3 = struct with fields:
    q1: [2x1 sym]
    q2: [2x1 sym]
    q3: [2x1 sym]
    q4: [2x1 sym]

```

```matlab
singular_4 = solve(det_4 == 0, [q1 q2 q3 q4])
```

```matlabTextOutput
singular_4 = struct with fields:
    q1: 0
    q2: 0
    q3: 0
    q4: 0

```

```matlab
j_s = subs(j, {q1,q2,q3,q4}, {0, pi/2, 0, 0})
```
j_s = 
 $\displaystyle \left(\begin{array}{cccc} -l_2  & -l_2  & 0 & -l_2 \newline l_1  & 0 & 1 & 0\newline 1 & 1 & 0 & 1 \end{array}\right)$
 

```matlab
rank(j_s) 
```

```matlabTextOutput
ans = 2
```
# planar 2R Robot
```matlab
syms q1 q2 dq1 dq2 'real'

% length l
l1=1
```

```matlabTextOutput
l1 = 1
```

```matlab
l2=1
```

```matlabTextOutput
l2 = 1
```

```matlab

% p in planar 2R
r=[l1*cos(q1)+l2*cos(q1+q2);
    l1*sin(q1)+l2*sin(q1+q2)]
```
r = 
 $\displaystyle \left(\begin{array}{c} \cos \left(q_1 +q_2 \right)+\cos \left(q_1 \right)\newline \sin \left(q_1 +q_2 \right)+\sin \left(q_1 \right) \end{array}\right)$
 

```matlab

% jacobian
J=jacobian(r,[q1 q2])
```
J = 
 $\displaystyle \left(\begin{array}{cc} -\sin \left(q_1 +q_2 \right)-\sin \left(q_1 \right) & -\sin \left(q_1 +q_2 \right)\newline \cos \left(q_1 +q_2 \right)+\cos \left(q_1 \right) & \cos \left(q_1 +q_2 \right) \end{array}\right)$
 

```matlab

% time derivative of the Jacobian
Jder=[-l1*cos(q1)*dq1-l2*cos(q1+q2)*(dq1+dq2) -l2*cos(q1+q2)*(dq1+dq2);
    -l1*sin(q1)*dq1-l2*sin(q1+q2)*(dq1+dq2) -l2*sin(q1+q2)*(dq1+dq2)];
```

```matlab
qa = [3*pi/4 -pi/2]'
```

```matlabTextOutput
qa = 2x1
    2.3562
   -1.5708

```

```matlab
qb = [pi/2 -pi/2]'
```

```matlabTextOutput
qb = 2x1
1.5708
   -1.5708

```

```matlab
FA=(10/sqrt(2))*[1 1]'
```

```matlabTextOutput
FA = 2x1
    7.0711
    7.0711

```

```matlab
j_qa = subs(J, [q1 q2], qa')
```
j_qa = 
 $\displaystyle \left(\begin{array}{cc} -\sqrt{2} & -\frac{\sqrt{2}}{2}\newline 0 & \frac{\sqrt{2}}{2} \end{array}\right)$
 

```matlab
tau_a = j_qa' * FA 
```
tau_a = 
 $\displaystyle \left(\begin{array}{c} -10\newline 0 \end{array}\right)$
 

```matlab

j_qb = subs(J, [q1 q2], qb')
```
j_qb = 
 $\displaystyle \left(\begin{array}{cc} -1 & 0\newline 1 & 1 \end{array}\right)$
 

```matlab
FB = (10/sqrt(2))*[-1 -1]'
```

```matlabTextOutput
FB = 2x1
   -7.0711
   -7.0711

```

```matlab
tau_b = -j_qb' * FB
```
tau_b = 
 $\displaystyle \left(\begin{array}{c} 0\newline 5\,\sqrt{2} \end{array}\right)$
 


```matlab
syms t T real
q1 = pi/4 + (pi/4)*(3*(t/T)^2-2*(t/T)^3)
```
q1 = 
 $\displaystyle \frac{\pi }{4}+\frac{\pi \,{\left(\frac{3\,t^2 }{T^2 }-\frac{2\,t^3 }{T^3 }\right)}}{4}$
 

```matlab
q1_dot = diff(q1, {t})
```
q1_dot = 
 $\displaystyle \frac{\pi \,{\left(\frac{6\,t}{T^2 }-\frac{6\,t^2 }{T^3 }\right)}}{4}$
 

```matlab
q1_dot_dot = diff(q1_dot, {t})
```
q1_dot_dot = 
 $\displaystyle -\frac{\pi \,{\left(\frac{12\,t}{T^3 }-\frac{6}{T^2 }\right)}}{4}$
 

```matlab

q2 = -pi/2*(1-cos(pi*t/T))
```
q2 = 
 $\displaystyle \frac{\pi \,{\left(\cos \left(\frac{\pi \,t}{T}\right)-1\right)}}{2}$
 

```matlab
q2_dot = diff(q2, {t})
```
q2_dot = 
 $\displaystyle -\frac{\pi^2 \,\sin \left(\frac{\pi \,t}{T}\right)}{2\,T}$
 

```matlab
q2_dot_dot = diff(q2_dot, {t})
```
q2_dot_dot = 
 $\displaystyle -\frac{\pi^3 \,\cos \left(\frac{\pi \,t}{T}\right)}{2\,T^2 }$
 

```matlab
q1_0 = subs(q1, {t}, {0})
```
q1_0 = 
 $\displaystyle \frac{\pi }{4}$
 

```matlab
q1_T = subs(q1, {t}, {T})
```
q1_T = 
 $\displaystyle \frac{\pi }{2}$
 

```matlab

q2_0 = subs(q2, {t}, {0})
```
q2_0 = 
 $\displaystyle 0$
 

```matlab
q2_T = subs(q2, {t}, {T})
```
q2_T = 
 $\displaystyle -\pi $
 

```matlab
q1_dot_0 = subs(q1_dot, {t}, {0})
```
q1_dot_0 = 
 $\displaystyle 0$
 

```matlab
q1_dot_T = subs(q1_dot, {t}, {T})
```
q1_dot_T = 
 $\displaystyle 0$
 

```matlab

q2_dot_0 = subs(q2_dot, {t}, {0})
```
q2_dot_0 = 
 $\displaystyle 0$
 

```matlab
q2_dot_T = subs(q2_dot, {t}, {T})
```
q2_dot_T = 
 $\displaystyle 0$
 

```matlab
q1_dot_dot_0 = subs(q1_dot_dot, {t}, {0})
```
q1_dot_dot_0 = 
 $\displaystyle \frac{3\,\pi }{2\,T^2 }$
 

```matlab
q1_dot_dot_T = subs(q1_dot_dot, {t}, {T})
```
q1_dot_dot_T = 
 $\displaystyle -\frac{3\,\pi }{2\,T^2 }$
 

```matlab

q2_dot_dot_0 = subs(q2_dot_dot, {t}, {0})
```
q2_dot_dot_0 = 
 $\displaystyle -\frac{\pi^3 }{2\,T^2 }$
 

```matlab
q2_dot_dot_T = subs(q2_dot_dot, {t}, {T})
```
q2_dot_dot_T = 
 $\displaystyle \frac{\pi^3 }{2\,T^2 }$
 

```matlab
q1_dot_dot_dot = diff(q1_dot_dot, {t})
```
q1_dot_dot_dot = 
 $\displaystyle -\frac{3\,\pi }{T^3 }$
 

```matlab
q2_dot_dot_dot = diff(q2_dot_dot, {t})
```
q2_dot_dot_dot = 
 $\displaystyle \frac{\pi^4 \,\sin \left(\frac{\pi \,t}{T}\right)}{2\,T^3 }$
 


```matlab
clf
T=5
```

```matlabTextOutput
T = 5
```

```matlab
hold on
fplot(subs(q1_dot_dot), [0,T])
hold off 
```

![figure_0.png](ex_2021mlx_media/figure_0.png)

```matlab

clf
T=5
```

```matlabTextOutput
T = 5
```

```matlab
hold on
fplot(subs(q2_dot_dot), [0,T])
hold off 
```

![figure_1.png](ex_2021mlx_media/figure_1.png)

```matlab
clf
T=5
```

```matlabTextOutput
T = 5
```

```matlab
hold on
fplot(subs(q1_dot), [0,T])
hold off 
```

![figure_2.png](ex_2021mlx_media/figure_2.png)

```matlab

clf
T=5
```

```matlabTextOutput
T = 5
```

```matlab
hold on
fplot(subs(q2_dot), [0,T])
hold off 
```

![figure_3.png](ex_2021mlx_media/figure_3.png)

```matlab
syms t T real
q1 = pi/4 + (pi/4)*(3*(t/T)^2-2*(t/T)^3)
```
q1 = 
 $\displaystyle \frac{\pi }{4}+\frac{\pi \,{\left(\frac{3\,t^2 }{T^2 }-\frac{2\,t^3 }{T^3 }\right)}}{4}$
 

```matlab
q1_dot = diff(q1, {t})
```
q1_dot = 
 $\displaystyle \frac{\pi \,{\left(\frac{6\,t}{T^2 }-\frac{6\,t^2 }{T^3 }\right)}}{4}$
 

```matlab
q1_dot_dot = diff(q1_dot, {t})
```
q1_dot_dot = 
 $\displaystyle -\frac{\pi \,{\left(\frac{12\,t}{T^3 }-\frac{6}{T^2 }\right)}}{4}$
 

```matlab

q2 = -pi/2*(1-cos(pi*t/T))
```
q2 = 
 $\displaystyle \frac{\pi \,{\left(\cos \left(\frac{\pi \,t}{T}\right)-1\right)}}{2}$
 

```matlab
q2_dot = diff(q2, {t})
```
q2_dot = 
 $\displaystyle -\frac{\pi^2 \,\sin \left(\frac{\pi \,t}{T}\right)}{2\,T}$
 

```matlab
q2_dot_dot = diff(q2_dot, {t})
```
q2_dot_dot = 
 $\displaystyle -\frac{\pi^3 \,\cos \left(\frac{\pi \,t}{T}\right)}{2\,T^2 }$
 

```matlab

subs(q1_dot, {t}, {T/2})
```
ans = 
 $\displaystyle \frac{3\,\pi }{8\,T}$
 

```matlab
subs(q2_dot, {t}, {T/2})
```
ans = 
 $\displaystyle -\frac{\pi^2 }{2\,T}$
 


```matlab
syms x gamma r real
eq = x^2 + (gamma*(x+0.8)^2)+1.1-r^2 == 0
```
eq = 
 $\displaystyle x^2 -r^2 +\gamma \,{{\left(x+\frac{4}{5}\right)}}^2 +\frac{11}{10}=0$
 

```matlab
sol = solve(eq, x, 'ReturnConditions',true)
```

```matlabTextOutput
sol = struct with fields:
             x: [2x1 sym]
    parameters: [1x0 sym]
    conditions: [2x1 sym]

```

```matlab
s = sol.x
```
s = 
 $\displaystyle \left(\begin{array}{c} -\frac{8\,\gamma -5\,\sqrt{4\,\gamma \,r^2 -\frac{174\,\gamma }{25}+4\,r^2 -\frac{22}{5}}}{10\,{\left(\gamma +1\right)}}\newline -\frac{8\,\gamma +5\,\sqrt{4\,\gamma \,r^2 -\frac{174\,\gamma }{25}+4\,r^2 -\frac{22}{5}}}{10\,{\left(\gamma +1\right)}} \end{array}\right)$
 

```matlab
sol.parameters
```

```matlabTextOutput
 
ans =
 
Empty sym: 1-by-0
 
```

```matlab
sol.conditions
```
ans = 
 $\displaystyle \left(\begin{array}{c} \gamma \not= -1\wedge 174\,\gamma +110\le r^2 \,{\left(100\,\gamma +100\right)}\newline \gamma \not= -1\wedge 174\,\gamma +110\le r^2 \,{\left(100\,\gamma +100\right)} \end{array}\right)$
 

```matlab
gamma = tan(deg2rad(-20))
```

```matlabTextOutput
gamma = -0.3640
```

```matlab
r = 0.5+0.4
```

```matlabTextOutput
r = 0.9000
```

```matlab
double(subs(s))
```

```matlabTextOutput
ans = 2x1
    0.8040
    0.1116

```

```matlab

x = double(subs(s(2)))
```

```matlabTextOutput
x = 0.1116
```

```matlab
pw = [x; double(subs(gamma*(x+0.8)+1.1))]
```

```matlabTextOutput
pw = 2x1
    0.1116
    0.7682

```

```matlab

prv_1 = pw(1)+0.6*sin(deg2rad(20))
```

```matlabTextOutput
prv_1 = 0.3168
```

```matlab
prv_2 = pw(2)-0.6*cos(deg2rad(20))
```

```matlabTextOutput
prv_2 = 0.2044
```

```matlab

l1=0.5
```

```matlabTextOutput
l1 = 0.5000
```

```matlab
l2=0.4
```

```matlabTextOutput
l2 = 0.4000
```

```matlab
%% Inverse kinematic
% x of p and y of p
px= double(prv_1)
```

```matlabTextOutput
px = 0.3168
```

```matlab
py= double(prv_2)
```

```matlabTextOutput
py = 0.2044
```

```matlab

%%  q2
% second joint computations
c2=(px^2+py^2-l1^2-l2^2)/(2*l1*l2);
    s2pos=sqrt(1-c2^2);%simplify(sqrt(1-c2^2));  
    %other solution: -sqrt(1-c2^2)
    s2neg= -sqrt(1-c2^2);%simplify(-sqrt(1-c2^2));  

%% q1
% first joint computations
detM=l1^2+l2^2+2*l1*l2*c2;

% positive solution of q1
s1pos=(py*(l1+l2*c2)-px*l2*s2pos)/detM; %simplify((py*(l1+l2*c2)-px*l2*s2pos)/detM);
c1pos=(px*(l1+l2*c2)+py*l2*s2pos)/detM;%simplify((px*(l1+l2*c2)+py*l2*s2pos)/detM);

% negative solution of q1
s1neg=(py*(l1+l2*c2)-px*l2*s2neg)/detM;%simplify((py*(l1+l2*c2)-px*l2*s2neg)/detM);
c1neg=(px*(l1+l2*c2)+py*l2*s2neg)/detM;%simplify((px*(l1+l2*c2)+py*l2*s2neg)/detM);

% output positive
fprintf("Positive configuration\n")
```

```matlabTextOutput
Positive configuration
```

```matlab
q01p=atan2(s1pos,c1pos);%simplify(atan2(s1pos,c1pos));
q02p=atan2(s2pos,c2);%simplify(atan2(s2pos,c2));
q0p=[q01p; q02p]
```

```matlabTextOutput
q0p = 2x1
   -0.3345
    2.3046

```

```matlab
    double(q0p)
```

```matlabTextOutput
ans = 2x1
   -0.3345
    2.3046

```

```matlab
%q0p=double([q01p; q02p])

fprintf("Negative configuration\n")
```

```matlabTextOutput
Negative configuration
```

```matlab
% output negative
q01n=atan2(s1neg,c1neg);%simplify(atan2(s1neg,c1neg));
q02n=atan2(s2neg,c2);%simplify(atan2(s2neg,c2));
q0n =[q01n;q02n]
```

```matlabTextOutput
q0n = 2x1
1.4805
   -2.3046

```

```matlab
    double(q0n)
```

```matlabTextOutput
ans = 2x1
1.4805
   -2.3046

```

```matlab
qrv = q0n
```

```matlabTextOutput
qrv = 2x1
1.4805
   -2.3046

```

```matlab

b = deg2rad(-20)
```

```matlabTextOutput
b = -0.3491
```

```matlab
v_vers = [cos(b); sin(b)]
```

```matlabTextOutput
v_vers = 2x1
    0.9397
   -0.3420

```

```matlab

v=v_vers*0.3
```

```matlabTextOutput
v = 2x1
    0.2819
   -0.1026

```

```matlab

syms q1 q2 dq1 dq2 'real'

% length l
l1=0.5
```

```matlabTextOutput
l1 = 0.5000
```

```matlab
l2=0.4
```

```matlabTextOutput
l2 = 0.4000
```

```matlab

% p in planar 2R
r=[l1*cos(q1)+l2*cos(q1+q2);
    l1*sin(q1)+l2*sin(q1+q2)]
```
r = 
 $\displaystyle \left(\begin{array}{c} \frac{2\,\cos \left(q_1 +q_2 \right)}{5}+\frac{\cos \left(q_1 \right)}{2}\newline \frac{2\,\sin \left(q_1 +q_2 \right)}{5}+\frac{\sin \left(q_1 \right)}{2} \end{array}\right)$
 

```matlab

% jacobian
J=jacobian(r,[q1 q2])
```
J = 
 $\displaystyle \left(\begin{array}{cc} -\frac{2\,\sin \left(q_1 +q_2 \right)}{5}-\frac{\sin \left(q_1 \right)}{2} & -\frac{2\,\sin \left(q_1 +q_2 \right)}{5}\newline \frac{2\,\cos \left(q_1 +q_2 \right)}{5}+\frac{\cos \left(q_1 \right)}{2} & \frac{2\,\cos \left(q_1 +q_2 \right)}{5} \end{array}\right)$
 

```matlab

q1 = qrv(1)
```

```matlabTextOutput
q1 = 1.4805
```

```matlab
q2 = qrv(2)
```

```matlabTextOutput
q2 = -2.3046
```

```matlab
j_rv = subs(J)
```
j_rv = 
 $\displaystyle \left(\begin{array}{cc} \frac{2\,\sin \left(\frac{1855628870113505}{2251799813685248}\right)}{5}-\frac{\sin \left(\frac{3333766676176955}{2251799813685248}\right)}{2} & \frac{2\,\sin \left(\frac{1855628870113505}{2251799813685248}\right)}{5}\newline \frac{\cos \left(\frac{3333766676176955}{2251799813685248}\right)}{2}+\frac{2\,\cos \left(\frac{1855628870113505}{2251799813685248}\right)}{5} & \frac{2\,\cos \left(\frac{1855628870113505}{2251799813685248}\right)}{5} \end{array}\right)$
 

```matlab
double(j_rv)
```

```matlabTextOutput
ans = 2x2
   -0.2044    0.2936
    0.3168    0.2717

```

```matlab

A = j_rv
```
A = 
 $\displaystyle \left(\begin{array}{cc} \frac{2\,\sin \left(\frac{1855628870113505}{2251799813685248}\right)}{5}-\frac{\sin \left(\frac{3333766676176955}{2251799813685248}\right)}{2} & \frac{2\,\sin \left(\frac{1855628870113505}{2251799813685248}\right)}{5}\newline \frac{\cos \left(\frac{3333766676176955}{2251799813685248}\right)}{2}+\frac{2\,\cos \left(\frac{1855628870113505}{2251799813685248}\right)}{5} & \frac{2\,\cos \left(\frac{1855628870113505}{2251799813685248}\right)}{5} \end{array}\right)$
 

```matlab
B = v
```

```matlabTextOutput
B = 2x1
    0.2819
   -0.1026

```

```matlab
q_rv_dot = linsolve(A,B)
```
q_rv_dot = 
 $\displaystyle \begin{array}{l} \left(\begin{array}{c} -\frac{20313596816708264\,\cos \left(\frac{1855628870113505}{2251799813685248}\right)+7393544592166489\,\sin \left(\frac{1855628870113505}{2251799813685248}\right)}{36028797018963968\,\sigma_1 }\newline \frac{101567984083541320\,\cos \left(\frac{3333766676176955}{2251799813685248}\right)-36967722960832445\,\sin \left(\frac{3333766676176955}{2251799813685248}\right)+81254387266833056\,\cos \left(\frac{1855628870113505}{2251799813685248}\right)+29574178368665956\,\sin \left(\frac{1855628870113505}{2251799813685248}\right)}{144115188075855872\,\sigma_1 } \end{array}\right)\\\mathrm{}\\\textrm{where}\\\mathrm{}\\\;\;\sigma_1 =\cos \left(\frac{3333766676176955}{2251799813685248}\right)\,\sin \left(\frac{1855628870113505}{2251799813685248}\right)+\sin \left(\frac{3333766676176955}{2251799813685248}\right)\,\cos \left(\frac{1855628870113505}{2251799813685248}\right)\end{array}$
 

```matlab
double(q_rv_dot)
```

```matlabTextOutput
ans = 2x1
   -0.7185
    0.4601

```

```matlab
% cubic polynomial
```
# Q1
```matlab
syms tau real
syms a b c d real
syms Dq real
syms qin qfin real
syms vin vfin real
syms T real

qn = a*tau^3 + b*tau^2 +c*tau + d 
```
qn = 
 $\displaystyle a\,\tau^3 +b\,\tau^2 +c\,\tau +d$
 

```matlab
q = qin + Dq*qn
```
q = 
 $\displaystyle \textrm{qin}+\textrm{Dq}\,{\left(a\,\tau^3 +b\,\tau^2 +c\,\tau +d\right)}$
 

```matlab

eq_1 = a+b+c==1
```
eq_1 = 
 $\displaystyle a+b+c=1$
 

```matlab
eq_2 = c==vin*T/Dq 
```
eq_2 = 
 $\displaystyle c=\frac{T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab
eq_3 = 3*a+2*b+c == vfin*T/Dq
```
eq_3 = 
 $\displaystyle 3\,a+2\,b+c=\frac{T\,\textrm{vfin}}{\textrm{Dq}}$
 

```matlab

s = solve([eq_1, eq_2, eq_3], [a, b, c])
```

```matlabTextOutput
s = struct with fields:
    a: (T*vfin)/Dq + (T*vin)/Dq - 2
    b: 3 - (2*T*vin)/Dq - (T*vfin)/Dq
    c: (T*vin)/Dq

```

```matlab
a = simplify(s.a)
```
a = 
 $\displaystyle \frac{T\,\textrm{vfin}-2\,\textrm{Dq}+T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab
b = simplify(s.b)
```
b = 
 $\displaystyle -\frac{T\,\textrm{vfin}-3\,\textrm{Dq}+2\,T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab
c = simplify(s.c)
```
c = 
 $\displaystyle \frac{T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab

%input
T = 2
```

```matlabTextOutput
T = 2
```

```matlab

qin = pi
```

```matlabTextOutput
qin = 3.1416
```

```matlab
qfin = 1.4805
```

```matlabTextOutput
qfin = 1.4805
```

```matlab

vin = 0
```

```matlabTextOutput
vin = 0
```

```matlab
vfin = -0.7185
```

```matlabTextOutput
vfin = -0.7185
```

```matlab

Dq = double(qfin-qin);
a = double(subs(a))
```

```matlabTextOutput
a = -1.1349
```

```matlab
b = double(subs(b))
```

```matlabTextOutput
b = 2.1349
```

```matlab
c = double(subs(c))
```

```matlabTextOutput
c = 0
```

```matlab
d = 0
```

```matlabTextOutput
d = 0
```

```matlab

syms t real
q_t = subs(q, {tau}, {t/T})
```
q_t = 
 $\displaystyle \textrm{qin}+\textrm{Dq}\,{\left(\frac{a\,t^3 }{8}+\frac{b\,t^2 }{4}+\frac{c\,t}{2}+d\right)}$
 

```matlab
q_t = simplify(subs(q_t));

clf
hold on
fplot(q_t, [0,T])
hold off 
```

![figure_4.png](ex_2021mlx_media/figure_4.png)

```matlab
q_dot = diff(q_t, {t})
```
q_dot = 
 $\displaystyle \frac{14338537717125850958202568576401\,t^2 }{20282409603651670423947251286016}-\frac{8990882771090890996441888283099\,t}{5070602400912917605986812821504}$
 

```matlab

clf
hold on
fplot(q_dot, [0,T])
hold off 
```

![figure_5.png](ex_2021mlx_media/figure_5.png)
# Q2
```matlab
% cubic polynomial
syms tau real
syms a b c d real
syms Dq real
syms qin qfin real
syms vin vfin real
syms T real

qn = a*tau^3 + b*tau^2 +c*tau + d 
```
qn = 
 $\displaystyle a\,\tau^3 +b\,\tau^2 +c\,\tau +d$
 

```matlab
q = qin + Dq*qn
```
q = 
 $\displaystyle \textrm{qin}+\textrm{Dq}\,{\left(a\,\tau^3 +b\,\tau^2 +c\,\tau +d\right)}$
 

```matlab

eq_1 = a+b+c==1
```
eq_1 = 
 $\displaystyle a+b+c=1$
 

```matlab
eq_2 = c==vin*T/Dq 
```
eq_2 = 
 $\displaystyle c=\frac{T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab
eq_3 = 3*a+2*b+c == vfin*T/Dq
```
eq_3 = 
 $\displaystyle 3\,a+2\,b+c=\frac{T\,\textrm{vfin}}{\textrm{Dq}}$
 

```matlab

s = solve([eq_1, eq_2, eq_3], [a, b, c])
```

```matlabTextOutput
s = struct with fields:
    a: (T*vfin)/Dq + (T*vin)/Dq - 2
    b: 3 - (2*T*vin)/Dq - (T*vfin)/Dq
    c: (T*vin)/Dq

```

```matlab
a = simplify(s.a)
```
a = 
 $\displaystyle \frac{T\,\textrm{vfin}-2\,\textrm{Dq}+T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab
b = simplify(s.b)
```
b = 
 $\displaystyle -\frac{T\,\textrm{vfin}-3\,\textrm{Dq}+2\,T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab
c = simplify(s.c)
```
c = 
 $\displaystyle \frac{T\,\textrm{vin}}{\textrm{Dq}}$
 

```matlab

%input
T = 2
```

```matlabTextOutput
T = 2
```

```matlab

qin = 0
```

```matlabTextOutput
qin = 0
```

```matlab
qfin = -2.3016
```

```matlabTextOutput
qfin = -2.3016
```

```matlab

vin = 0
```

```matlabTextOutput
vin = 0
```

```matlab
vfin = 0.4601
```

```matlabTextOutput
vfin = 0.4601
```

```matlab

Dq = double(qfin-qin);
a = double(subs(a))
```

```matlabTextOutput
a = -2.3998
```

```matlab
b = double(subs(b))
```

```matlabTextOutput
b = 3.3998
```

```matlab
c = double(subs(c))
```

```matlabTextOutput
c = 0
```

```matlab
d = 0
```

```matlabTextOutput
d = 0
```

```matlab

syms t real
q_t = subs(q, {tau}, {t/T})
```
q_t = 
 $\displaystyle \textrm{qin}+\textrm{Dq}\,{\left(\frac{a\,t^3 }{8}+\frac{b\,t^2 }{4}+\frac{c\,t}{2}+d\right)}$
 

```matlab
q_t = simplify(subs(q_t));

clf
hold on
fplot(q_t, [0,T])
hold off 
```

![figure_6.png](ex_2021mlx_media/figure_6.png)

```matlab
q_dot = diff(q_t, {t})
```
q_dot = 
 $\displaystyle \frac{t\,{\left(27617\,t-78250\right)}}{20000}+\frac{27617\,t^2 }{40000}$
 

```matlab

clf
hold on
fplot(q_dot, [0,T])
hold off 
```

![figure_7.png](ex_2021mlx_media/figure_7.png)
