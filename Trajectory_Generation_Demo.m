clc
clear variables

% Changable Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
poc = [165,865,-34.2];  % X-Y-Z position of ceneter 
stepin = 0.2;  % Step-in depth
arcpn = 32;  % Arc Point Number
linpd = 2;  % Distances betweem Linear Points
l=180;  % Length of side
d=40;  % Depth of pocket
s=37;  % Shrink of pocket
r=15;  % Radius of corner
tool_length=305.5;  % Tool Length
tooloffset=16;  % Offset of the tool (Positive for longer)
tool_radius=7.5;  % Radius of the tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loopnum=ceil(d/stepin);

% Store point data
PD=zeros(arcpn*4*(2+loopnum),3);
PC=0;  % Point counter

% Calculate dTCP: horizontal distance between TCP to the workpiece
alpha=atan(d/s);
dTCP=tool_radius*(1-cos(alpha))/sin(alpha);

% Record Belongings of points in different portions, 0 for upper and lower
SidePointer=ones(71024,1);

%% First Loop
zpos=poc(3);
for arc=1:4
    % Arc Portion
    if arc==1
        arcpos=[1,1];
    elseif arc==2
        arcpos=[-1,1];
    elseif arc==3
        arcpos=[-1,-1];
    elseif arc==4
        arcpos=[1,-1];
    end    
    for oln=0:arcpn-1
        xypos=[poc(1),poc(2)]+(l-2*r)*(0.5)*arcpos+(r-dTCP)*[cos(0.5*pi*((oln)/(arcpn-1)+arc-1)),sin(0.5*pi*((oln)/(arcpn-1)+arc-1))];
        
        % Store data
        PC=PC+1;
        PD(PC,:)=[xypos,zpos];
    end
    
    % Linear Portion 
    linpn = (l-2*r)/linpd;  % Calculate linear point number
    for oln=1:linpn-1
        if arc==1
            xypos=[poc(1),poc(2)]+[(l-2*r)*(0.5)-oln*linpd,l*(0.5)-dTCP];
        elseif arc==2
            xypos=[poc(1),poc(2)]+[-l*(0.5)+dTCP,(l-2*r)*(0.5)-oln*linpd];
        elseif arc==3
            xypos=[poc(1),poc(2)]+[-(l-2*r)*(0.5)+oln*linpd,-l*(0.5)+dTCP];
        elseif arc==4
            xypos=[poc(1),poc(2)]+[l*(0.5)-dTCP,-(l-2*r)*(0.5)+oln*linpd];
        end

        % Store data
        PC=PC+1;
        PD(PC,:)=[xypos,zpos];
        if arc==2 || arc==4
            SidePointer(PC,1)=0;
        else
            SidePointer(PC,1)=2;            
        end
    end
end

%% Middle Loops
for loop=1:loopnum
    tanthz=(stepin)/(4*(l-2*(s/d)*stepin*loop-2*r)+2*pi*r);  % tan(th)
    tanthxy=(stepin*(s/d))/(4*(l-2*(s/d)*stepin*loop-2*r));  % tan(th_xy)

    for arc=1:4
        % Arc Portion
        if arc==1
            arcpos=[1,1,0];
        elseif arc==2
            arcpos=[-1,1,0];
        elseif arc==3
            arcpos=[-1,-1,0];
        elseif arc==4
            arcpos=[1,-1,0];
        end
        for oln=0:arcpn-1
            poc2arccenter_xy = ((l-2*r)*(0.5)-(0.25)*(4*(loop-1)+(arc-1))*stepin*(s/d))*arcpos;
            arccenter2arc_xy = (r-dTCP)*[cos(0.5*pi*((oln)/(arcpn-1)+arc-1)),sin(0.5*pi*((oln)/(arcpn-1)+arc-1)),0];
            depthin = ((0.25)*(4*(loop-1)+(arc-1))*stepin + 0.5*pi*r*(oln/(arcpn-1))*tanthz)*[0,0,1];
            xyzpos = poc + poc2arccenter_xy + arccenter2arc_xy - depthin;
            
            % Store data
            PC=PC+1;
            PD(PC,:)=xyzpos;
        end

        % Linear Portion
        linpn=ceil((l-2*(s/d)*stepin*loop-2*r)/linpd);
        for oln=1:linpn-1
            if arc==1  % Right side
                poc2line_xy=[poc(1),poc(2),0]+[(l-2*(s/d)*stepin*loop-2*r)*(0.5)-oln*linpd,(l-2*(s/d)*stepin*loop)*(0.5)-dTCP-oln*linpd*tanthxy-0.25*(arc-5)*stepin*(s/d),0];
                depthin= ((0.25)*(4*(loop-1)+(arc-1))*stepin + (0.5*pi*r+oln*linpd)*tanthz)*[0,0,1];
            elseif arc==2  % Up side
                poc2line_xy=[poc(1),poc(2),0]+[-(l-2*(s/d)*stepin*loop)*(0.5)+dTCP+oln*linpd*tanthxy+0.25*(arc-5)*stepin*(s/d),(l-2*(s/d)*stepin*loop-2*r)*(0.5)-oln*linpd,0];
                depthin= ((0.25)*(4*(loop-1)+(arc-1))*stepin + (0.5*pi*r+oln*linpd)*tanthz)*[0,0,1];
            elseif arc==3  % Left side
                poc2line_xy=[poc(1),poc(2),0]+[-(l-2*(s/d)*stepin*loop-2*r)*(0.5)+oln*linpd,-(l-2*(s/d)*stepin*loop)*(0.5)+dTCP+oln*linpd*tanthxy+0.25*(arc-5)*stepin*(s/d),0];
                depthin= ((0.25)*(4*(loop-1)+(arc-1))*stepin + (0.5*pi*r+oln*linpd)*tanthz)*[0,0,1];
            elseif arc==4  % Down side
                poc2line_xy=[poc(1),poc(2),0]+[(l-2*(s/d)*stepin*loop)*(0.5)-dTCP-oln*linpd*tanthxy-0.25*(arc-5)*stepin*(s/d),-(l-2*(s/d)*stepin*loop-2*r)*(0.5)+oln*linpd,0];
                depthin= ((0.25)*(4*(loop-1)+(arc-1))*stepin + (0.5*pi*r+oln*linpd)*tanthz)*[0,0,1];
            end
            xyzpos=poc2line_xy+[0,0,poc(3)]-depthin;

            % Store data
            PC=PC+1;
            PD(PC,:)=xyzpos;
            if arc==2 || arc==4
                SidePointer(PC,1)=0;
            else
                SidePointer(PC,1)=2;                        
            end
        end
    end
end

%% Last Loop
zpos=poc(3)-40;
for arc=1:4
    % Arc Portion
    if arc==1
        arcpos=[1,1];
    elseif arc==2
        arcpos=[-1,1];
    elseif arc==3
        arcpos=[-1,-1];
    elseif arc==4
        arcpos=[1,-1];
    end
    for oln=0:arcpn-1
        xypos=[poc(1),poc(2)]+(l-2*r-2*s)*(0.5)*arcpos+(r-dTCP)*[cos(0.5*pi*((oln)/(arcpn-1)+arc-1)),sin(0.5*pi*((oln)/(arcpn-1)+arc-1))];
        % Store data
        PC=PC+1;
        PD(PC,:)=[xypos,zpos];
    end

    % Linear Portion 
    linpn = (l-2*r-2*s)/linpd;  % Calculate linear point number
    for oln=1:linpn-1
        if arc==1
            xypos=[poc(1),poc(2)]+[(l-2*r-2*s)*(0.5)-oln*linpd,(l-2*s)*(0.5)-dTCP];
        elseif arc==2
            xypos=[poc(1),poc(2)]+[-(l-2*s)*(0.5)+dTCP,(l-2*r-2*s)*(0.5)-oln*linpd];
        elseif arc==3
            xypos=[poc(1),poc(2)]+[-(l-2*r-2*s)*(0.5)+oln*linpd,-(l-2*s)*(0.5)+dTCP];
        elseif arc==4
            xypos=[poc(1),poc(2)]+[(l-2*s)*(0.5)-dTCP,-(l-2*r-2*s)*(0.5)+oln*linpd];
        end

        % Store data
        PC=PC+1;
        PD(PC,:)=[xypos,zpos];
        if arc==2 || arc==4
            SidePointer(PC,1)=0;
        else
            SidePointer(PC,1)=2;  
        end
    end
end
PD(:,3)=PD(:,3)-tooloffset;

%% Transform to Joint Space
% Global frame to tool frame
Rt2g=[0  -1   0;
      0   0   1;
     -1   0   0];
Pt2g=[871.6 -1850 1914.4];

% Transform Point Data to Global frame
PDG=zeros(PC,3);
for i=1:PC
    PDG(i,:)=Pt2g+transpose(Rt2g*transpose(PD(i,:)));
end

% Create Joint Data based on Inverse kinematics
PDJ=zeros(PC,6);
t=[1, 0, 0];
b=[0, 0, 1];
n=[0, -1, 0];

for i=1:PC
    PDJ(i,:)=(360/(2*pi))*inv_kine(n,t,b,PDG(i,:),tool_length);
    if PDJ(i,6)>0
        PDJ(i,4)=PDJ(i,4)-180;
        PDJ(i,5)=-PDJ(i,5);
        PDJ(i,6)=PDJ(i,6)-180;
    end
end

%% Show points in 3D plot
figure
%plot3(PD(:,1),PD(:,2),PD(:,3),'black');
hold on
for i=1:71024
    if SidePointer(i,1)==1
        plot3(PD(i,1),PD(i,2),PD(i,3),'black.');
    elseif SidePointer(i,1)==0
        plot3(PD(i,1),PD(i,2),PD(i,3),'red.');
    else
        plot3(PD(i,1),PD(i,2),PD(i,3),'blue.');
    end
end
xlabel('x')
xlim([0 300])
ylabel('y')
ylim([700 1000])
zlabel('z')
zlim([-100 -40])
hold off


%% Prepare MOD Files
% File 1
fid=fopen('SpeedVariation_1.mod','wt');

fprintf(fid,'%%%%%%\n');
fprintf(fid,'  VERSION:1\n');
fprintf(fid,'LANGUAGE:ENGLISH\n');
fprintf(fid,'%%%%%%\n\n');
fprintf(fid,'MODULE SpeedVariation_1\n');
fprintf(fid,'  PROC Proc_SpeedVariation_1()\n');
fprintf(fid,'    ! Start position = (-72.8553,-14.4801,27.1619,54.5626,-21.2122,-32.6507,9E9,9E9,9E9,9E9,9E9,9E9)\n    ! Tool Number    = 1\n    ! Spindle Speed  = 0 RPM\n');
fprintf(fid,'    MoveJ [[%f,%f,150.000],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v30,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(1,1),PD(1,2));
fprintf(fid,'    ! First Toolpath Point\n');
fprintf(fid,'    MoveL [[%f,%f,150.000],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v30,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(1,1),PD(1,2));
fprintf(fid,'    MoveL [[%f,%f,100.000],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v30,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(1,1),PD(1,2));
fprintf(fid,'    ! Plunge Move Starts\n');
fprintf(fid,'    MoveL [[%f,%f,%f],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v30,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(1,1),PD(1,2),(PD(1,3)+100.2));
fprintf(fid,'    ! Cutting Move Starts\n');
for i=1:49000
    if SidePointer(i,1)==0
        fprintf(fid,'    MoveL [[%f,%f,%f],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v10,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(i,1),PD(i,2),PD(i,3));
    elseif SidePointer(i,1)==1
        fprintf(fid,'    MoveL [[%f,%f,%f],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v30,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(i,1),PD(i,2),PD(i,3));
    else
        fprintf(fid,'    MoveL [[%f,%f,%f],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v40,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(i,1),PD(i,2),PD(i,3));
    end
end
fprintf(fid,'  ENDPROC\n');
fprintf(fid,'ENDMODULE\n');

% File 2
fid=fopen('SpeedVariation_2.mod','wt');
fprintf(fid,'%%%%%%\n');
fprintf(fid,'  VERSION:1\n');
fprintf(fid,'LANGUAGE:ENGLISH\n');
fprintf(fid,'%%%%%%\n\n');
fprintf(fid,'MODULE SpeedVariation_2\n');
fprintf(fid,'  PROC Proc_SpeedVariation_2()\n');
for i=49001:PC
    if SidePointer(i,1)==0
        fprintf(fid,'    MoveL [[%f,%f,%f],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v10,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(i,1),PD(i,2),PD(i,3));
    elseif SidePointer(i,1)==1
        fprintf(fid,'    MoveL [[%f,%f,%f],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v30,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(i,1),PD(i,2),PD(i,3));
    else
        fprintf(fid,'    MoveL [[%f,%f,%f],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v40,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(i,1),PD(i,2),PD(i,3));
    end
end

fprintf(fid,'    ! Cutting Move Ends\n');
fprintf(fid,'    MoveL [[%f,%f,150.000],[0.00030000,0.00000000,0.00000000,0.99999996],[-1,0,-1,1],[9E9,9E9,9E9,9E9,9E9,9E9]],v30,z1,tAutodesk1\\WObj:=wAutodesk1;\n',PD(i,1),PD(i,2));
fprintf(fid,'  ENDPROC\n');
fprintf(fid,'ENDMODULE\n');
fclose(fid);

%% Function Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matrix=rot_dh(al, th)
    % Calculate rotation matrix with DH-parameters: alpha, theta.
    matrix = [cos(th), -sin(th)*cos(al),  sin(th)*sin(al);
              sin(th),  cos(th)*cos(al), -cos(th)*sin(al);
              0,        sin(al),          cos(al)];
end

function j = inv_kine(n,t,b,p,tool_length)
    % DH_parameter and else for IK %%%%%%%%%%%%%%%%%%%%
    DH_a = [410 1075 165 0 0 0];
    DH_alpha = [-pi/2 0 -pi/2 pi/2 -pi/2 0];
    DH_d = [780 0 0 1056 0 250+tool_length];
    DH_theta = [0 -pi/2 0 0 0 0];
    arm_forward = true;
    bend_up = true;
    rev456 = false;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v1 = p - (DH_d(5)+DH_d(6))*n;

    % Shoulder singularity when the wrist is just above J1.
    if v1(1) == 0 && v1(2) == 0
        disp('Warning: Shoulder singularity.')
        j1 = 0;  % Keep J1 at initial position
    elseif arm_forward==true
        j1 = atan2(v1(2), v1(1));
    else
        j1 = atan2(-v1(2), -v1(1));
    end

    % Calculate J2 and J3 with v2 point from center of J2 to J5
    v2 = [0, 0, 0];
    if arm_forward==true
        v2(1) = v1(1)*cos(j1) + v1(2)*sin(j1) - DH_a(1);  % coordinate on x-y plane
    else
        v2(1) = -v1(1)*cos(j1) - v1(2) * sin(j1) - DH_a(1);
    end
    v2(3) = v1(3) - DH_d(1);  % coordinate on z plane
    
    dv2 = sqrt(v2(1)*v2(1) + v2(3)*v2(3));  % Length of v2
    l23 = DH_a(2);  % Length between J2 and J3
    l35 = sqrt(DH_a(3)*DH_a(3) + DH_d(4)*DH_d(4));  % Length between J3 and J5

    if dv2 >= (l23+l35)
        disp("Warning: elbow singularity.")
    end
   
    angle_325 = acos((l23*l23 + dv2*dv2 - l35*l35) / (2*l23*dv2));  % Angle J3-J2-J5
    angle_235 = acos((l23*l23 + l35*l35 - dv2*dv2) / (2*l23*l35));  % Angle J2-J3-J5
    angle_52h = atan2(v2(3), v2(1));  % Angle J5-J2-horizontal
    offset_35 = atan2(DH_a(3), DH_d(4));  % Offset of J3-J5 due to a(3)

    % Determine J2 and J3 based on specified configuration
    if arm_forward==true
        if bend_up==true
            j2 = pi/2 - angle_52h - angle_325;
            j3 = pi/2 + offset_35 - angle_235;
        else
            j2 = pi/2 - angle_52h + angle_325;
            j3 = -pi*(3/2) + offset_35 + angle_235;
        end
    else
        if bend_up==true
            j2 = pi/2 + angle_52h + angle_325;
            j3 = -pi*(3/2) + offset_35 + angle_235;
        else
            j2 = pi/2 + angle_52h - angle_325;
            j3 = -pi/2 - offset_35 + angle_235;
        end
    end

    % Calculate J4, J5 and J6
    R01 = rot_dh(DH_alpha(1), j1+DH_theta(1));
    R12 = rot_dh(DH_alpha(2), j2+DH_theta(2));
    R23 = rot_dh(DH_alpha(3), j3+DH_theta(3));
    R03 = R01*R12*R23;
    R06=[transpose(t) transpose(b) transpose(n)];
    
    R36=R03\R06;
    
    c5 = R36(3,3);  % cos(theta5)
    s4s6 = (R36(1,1) - c5 * R36(2,2)) / (c5*c5 - 1);  % sin(theta4)*sin(theta6)
    c4c6 = (c5 * R36(1,1) - R36(2,2)) / (c5*c5 - 1);  % cos(theta4)*cos(theta6)
    c4s6 = (R36(2,1) + c5 * R36(1,2)) / (-c5*c5 + 1);  % cos(theta4)*sin(theta6)
    s4c6 = (c5 * R36(2,1) + R36(1,2)) / (c5*c5 - 1);  % sin(theta4)*cos(theta6)


    s4plus6 = c4s6 + s4c6;  % sin(theta4+theta6)
    s4minus6 = s4c6 - c4s6;  % sin(theta4-theta6)
    c4plus6 = c4c6 - s4s6;  % cos(theta4+theta6)
    c4minus6 = c4c6 + s4s6;  % cos(theta4-theta6)

    ang4plus6 = atan2(s4plus6, c4plus6);  % theta4 + theta6
    ang4minus6 = atan2(s4minus6, c4minus6);  % theta4 - theta6

    j4 = 0.5*(ang4plus6 + ang4minus6);  % theta 4
    j6 = 0.5*(ang4plus6 - ang4minus6);  % theta 6

    s5 = R36(3,1)/cos(j6);  % sin(theta5)
    j5 = atan2(s5, c5);  % theta5

    if j4<0
        rev456=true;
    end
    
    if rev456
        if j4>0
            j4 = j4-pi;
        else
            j4 = j4 + pi;
        end
        j5 = -j5;
        if j6 > 0
            j6 = j6-pi;
        else
            j6 = j6 + pi;
        end
    end

    j=[j1 j2 j3 j4 j5 j6];
end

