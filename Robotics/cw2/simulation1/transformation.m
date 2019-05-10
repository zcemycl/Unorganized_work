
function [r1T, r2T, r3T, r4T] = transformation(Phi, L, mode)
phi1_T = Phi(1); phi2_T = Phi(2); phi3_T = Phi(3);
A0T = [1,0,L;
       0,1,0;
       0,0,1];
A10T = [1 , 0 ,  0; ...
        0 , 1 , L; ...
        0 , 0 , 1 ];
if mode == 'r'
    A21T = [cos(phi1_T) , sin(phi1_T) , L*sin(phi1_T) ; ...
           -sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
            0        , 0          , 1  ];

    A32T = [cos(phi2_T) , sin(phi2_T) , L*sin(phi2_T); ...
           -sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
           0         , 0         , 1 ];
    A43T = [cos(phi3_T) , sin(phi3_T) , L*sin(phi3_T); ...
           -sin(phi3_T) ,  cos(phi3_T) , L*cos(phi3_T); ...
            0         , 0         , 1 ];
elseif mode == 'f'
    A21T = [cos(phi1_T) , -sin(phi1_T) , -L*sin(phi1_T) ; ...
            sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
            0        , 0          , 1  ];

    A32T = [cos(phi2_T) , -sin(phi2_T) , -L*sin(phi2_T); ...
            sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
            0         , 0         , 1 ];
    A43T = [cos(phi3_T) , -sin(phi3_T) , -L*sin(phi3_T); ...
            sin(phi3_T) ,  cos(phi3_T) , L*cos(phi3_T); ...
            0         , 0         , 1 ];
end
if mode == 'r'
    r1T = A10T(:,3);   %position link 1
    r1T = r1T(1:2);
    r2T = A10T*A21T;   %position link 2
    r2T = r2T(:,3);
    r2T = r2T(1:2);
    r3T = A10T*A21T*A32T;
    r3T = r3T(:,3);
    r3T = r3T(1:2);
    r4T = A10T*A21T*A32T*A43T;
    r4T = r4T(:,3);
    r4T = r4T(1:2);
elseif mode == 'f'
    r1T = A10T*A43T*A32T*A21T;   %position link 1
    r1T = r1T(:,3);
    r1T = r1T(1:2);
    r2T = A10T*A43T*A32T;   %position link 2
    r2T = r2T(:,3);
    r2T = r2T(1:2);
    r3T = A10T*A43T;
    r3T = r3T(:,3);
    r3T = r3T(1:2);
    r4T = A10T;
    r4T = r4T(:,3);
    r4T = r4T(1:2);
end
end
