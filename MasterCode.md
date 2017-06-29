# MAE154B
Structural analysis code for the MAE 154B - Design of Aircraft Structures course.

%% STRUCTURAL ANALYSIS MASTER CODE FEAT BILL, MOHAMMED AND MATT

%Acknowledgements: 
%Matt: Cars
%Mohammed: my family
%Bill: Friends and family

clearvars -except Monte mc Load iterations W FOS
close all
%% Preliminary geometry (construct airfoil, plot spar caps + stringers)
%clc; clear;
%stringer location 
c = 5; %ft
M = .02 * c;
P = .4;

x_front = 0:0.01:P*c;
x_aft = P*c+0.01:0.01:c;
x = horzcat(x_front, x_aft);

camber_front = (M/(P^2)) * (2 * P * (x_front/c) - (x_front/c).^2); %Construct front camber line
camber_aft = M/((1-P)^2) * (1 - 2 * P + 2 * P * (x_aft/c) - (x_aft/c).^2); %Construct back camber line
camber = horzcat(camber_front, camber_aft); %One camber line!
t = (3) * (0.2969 * (x/c).^ 0.5 - 0.126*(x/c) - 0.3516*(x/c).^2 + 0.2843*(x/c).^3 - 0.1015*(x/c).^4);

%Plot the NACA 2412 
plot(x, camber, 'k--', 'LineWidth', 0.25); %Sanity check!
hold on;
z_upper = camber + t;
z_lower = camber - t;
plot(x, z_upper, 'k', 'LineWidth', 2);
plot(x, z_lower, 'k', 'LineWidth', 2);

%% Plot the main spar caps and webs
%Need to concatenate y in case x for stringer 1 is 0.4*c

sc_1 = round(rand*1.50,2) + 0.50; %Spar Cap 1 at at this x-location; keeps this 0 < x <= 2
sc_2 = round(rand*1.75,2) + 2; %Spar Cap 2 at this x-location; keeps this 2 < x < 3.75

sc_1 = 1.15;
sc_2 = 2.55;

index_sc_1 = find(abs(x - sc_1) < 0.001); %Find the index at which x = sc_1
index_sc_2 = find(abs(x - sc_2) < 0.001); %Weird syntax bc of floating-pt error

sc_1zu = z_upper(index_sc_1); %Find the y corresponding to that x
sc_1zl = z_lower(index_sc_1);

sc_2zu = z_upper(index_sc_2); %Same as above but for str_2
sc_2zl = z_lower(index_sc_2);

sc = [sc_2 sc_2zu; sc_2 sc_2zl; sc_1 sc_1zl; sc_1 sc_1zu];

%Plot out the spar caps and the corresponding webs
plot(sc_1, sc_1zu, 'rs', 'LineWidth', 2); %Plotting 4 spar caps
plot(sc_1, sc_1zl, 'rs', 'LineWidth', 2);
plot(sc_2, sc_2zu, 'rs', 'LineWidth', 2);
plot(sc_2, sc_2zl, 'rs', 'LineWidth', 2);

plot([sc_1 sc_1], [sc_1zl sc_1zu], 'r'); %Plotting corresponding webs
plot([sc_2 sc_2], [sc_2zl sc_2zu], 'r');

%% Determine number of stringers and where they are placed
num_str_1 = randi(8); %Number of stringers on top-right web
num_str_2 = randi(8); %Number of stringers on bottom-right web
num_str_3 = randi(8); %Number of stringers on bottom-left
num_str_4 = randi(8); %Number of stringers on top-left

num_str_1 = 6;
num_str_2 = 6;
num_str_3 = 5;
num_str_4 = 5;

NS = num_str_1 + num_str_2 + num_str_3 + num_str_4 + 4; %Total number of stirngers

incr_1 = round((index_sc_2 - index_sc_1)/(num_str_1 + 1), 0); %Increment from index_1 to index_2 at which each stringer will be placed
incr_2 = round((index_sc_2 - index_sc_1)/(num_str_2 + 1), 0);
incr_3 = round(index_sc_1/(num_str_3 + 1), 0);
incr_4 = round(index_sc_1/(num_str_4 + 1), 0);

%LOCATIONS OF STRINGERS FOR SECTION 1
for i = 1:num_str_1
   str_1_x(i) = index_sc_1 + incr_1 * i; %Calculate out your increments
   str_1_z(i) = z_upper(str_1_x(i));
end

str_1_x = x(str_1_x);
str_1 = [transpose(str_1_x) transpose(str_1_z)]; %x and y locations of stringers in section 1

%LOCATIONS OF STRINGERS FOR SECTION 1
for i = 1:num_str_2
    str_2_x(i) = index_sc_1 + incr_2 * i;
    str_2_z(i) = z_lower(str_2_x(i));    
end

str_2_x = x(str_2_x);
str_2 = [transpose(str_2_x) transpose(str_2_z)]; %x and y locations of stringers in section 2

%LOCATIONS OF STRINGERS FOR SECTION 3
for i = 1:num_str_3
    str_3_x(i) = incr_3 * i;
    str_3_z(i) = z_lower(str_3_x(i));    
end

str_3_x = x(str_3_x);
str_3 = [transpose(str_3_x) transpose(str_3_z)]; %x and y locations of stringers in section 3

%LOCATIONS OF STRINGERS FOR SECTION 4
for i = 1:num_str_4
    str_4_x(i) = incr_4 * i;
    str_4_z(i) = z_upper(str_4_x(i));    
end

str_4_x = x(str_4_x);
str_4 = [transpose(str_4_x) transpose(str_4_z)]; %x and y locations of stringers in section 3

%Plot some stringers!
scatter(str_1(:,1), str_1(:,2), 'r', 'filled');
scatter(str_2(:,1), str_2(:,2), 'r', 'filled');
scatter(str_3(:,1), str_3(:,2), 'r', 'filled');
scatter(str_4(:,1), str_4(:,2), 'r', 'filled');
ylim([-2.5 2.5]);

% Section1 = struct('Num_Str', num_str_1 + 2, 'Areas', rand(num_str_1 + 2,1)/10, 'Loc_Str', vertcat([sc_1 sc_1zu], str_1, [sc_2 sc_2zu]), 'Unit_Areas', 'ft^2', 'Thickness', rand/10);
% Section2 = struct('Num_Str', num_str_2 + 2, 'Areas', rand(num_str_2 + 2,1)/10, 'Loc_Str', vertcat([sc_1 sc_1zl], str_2, [sc_2 sc_2zl]), 'Unit_Areas', 'ft^2', 'Thickness', rand/10);
% Section3 = struct('Num_Str', num_str_3, 'Areas', rand(num_str_3,1)/10, 'Loc_Str', str_3, 'Unit_Areas', 'ft^2', 'Thickness', rand); 
% Section4 = struct('Num_Str', num_str_4, 'Areas', rand(num_str_4,1)/10, 'Loc_Str', str_4, 'Unit_Areas', 'ft^2', 'Thickness', Section3.Thickness); 
% % SparCaps = struct - need to do this so we can assign a thickness to the
% SparCap = struct('Thickness1', rand, 'Thickness2', rand);

Section1 = struct('Num_Str', num_str_1 + 2, 'Areas', zeros(num_str_1 + 2,1) + 0.105, 'Loc_Str', vertcat([sc_1 sc_1zu], str_1, [sc_2 sc_2zu]), 'Unit_Areas', 'ft^2', 'Thickness', 0.025);
Section2 = struct('Num_Str', num_str_2 + 2, 'Areas', zeros(num_str_2 + 2,1) + 0.105, 'Loc_Str', vertcat([sc_1 sc_1zl], str_2, [sc_2 sc_2zl]), 'Unit_Areas', 'ft^2', 'Thickness', 0.025);
Section3 = struct('Num_Str', num_str_3, 'Areas', zeros(num_str_3,1) + 0.125, 'Loc_Str', str_3, 'Unit_Areas', 'ft^2', 'Thickness', 0.025); 
Section4 = struct('Num_Str', num_str_4, 'Areas', zeros(num_str_4,1) + 0.125, 'Loc_Str', str_4, 'Unit_Areas', 'ft^2', 'Thickness', Section3.Thickness); 
% SparCaps = struct - need to do this so we can assign a thickness to the
SparCap = struct('Thickness1', rand, 'Thickness2', rand);

%% Analysis - Find centroid of our airfoil and apply web idealization
sum_Areas = sum(Section1.Areas) + sum(Section2.Areas) + sum(Section3.Areas) + sum(Section4.Areas); %sum of stringer + sc areas

num_1_z = dot(Section1.Areas, Section1.Loc_Str(:,2));
num_2_z = dot(Section2.Areas, Section2.Loc_Str(:,2));
num_3_z = dot(Section3.Areas, Section3.Loc_Str(:,2));
num_4_z = dot(Section4.Areas, Section4.Loc_Str(:,2));

num_1_x = dot(Section1.Areas, Section1.Loc_Str(:,1));
num_2_x = dot(Section2.Areas, Section2.Loc_Str(:,1));
num_3_x = dot(Section3.Areas, Section3.Loc_Str(:,1));
num_4_x = dot(Section4.Areas, Section4.Loc_Str(:,1));

ZC = sum([num_1_z, num_2_z, num_3_z, num_4_z])/sum_Areas;
XC = sum([num_1_x, num_2_x, num_3_x, num_4_x])/sum_Areas;
scatter(XC, ZC, 'ko', 'filled');
% hold off;

%Apply the web thickness idealization
%Spar Caps
%SC_left = Area + B_below + str_left + str_right
%Section 4
%if num_str_4 == 1
%Rotely apply first approximation to first top stringer
%else
%for i = 2:how every many stringers in section 4
    %B1 = (i-1)
    %B2 = (i+1)
    %Section4.Stringer
%end


%Bottom Web
%Rotely apply approximation to first top stringer
%Rotely apply to the last?
%for 2 = number of stringers on bottom


%end

%Calculate inertia terms
dz_1 = Section1.Loc_Str(:,2) - ZC;
dz_2 = Section2.Loc_Str(:,2) - ZC;
dz_3 = Section3.Loc_Str(:,2) - ZC;
dz_4 = Section4.Loc_Str(:,2) - ZC;

dx_1 = Section1.Loc_Str(:,1) - XC;
dx_2 = Section2.Loc_Str(:,1) - XC;
dx_3 = Section3.Loc_Str(:,1) - XC;
dx_4 = Section4.Loc_Str(:,1) - XC;

Ix = dot(Section1.Areas, dz_1.^2) + dot(Section2.Areas, dz_2.^2) + dot(Section3.Areas, dz_3.^2) + dot(Section4.Areas, dz_4.^2);
Iz = dot(Section1.Areas, dx_1.^2) + dot(Section2.Areas, dx_2.^2) + dot(Section3.Areas, dx_3.^2) + dot(Section4.Areas, dx_4.^2);
Ixz = dot(Section1.Areas, dx_1 .* dz_1) + dot(Section2.Areas, dx_2 .* dz_2) + dot(Section3.Areas, dx_3 .* dz_3) + dot(Section4.Areas, dx_4 .* dz_4);

%% Analysis - Calculation of dP, shear center, and shear flows
for step = 1:3
    if step == 1
        Vx = 1; Vz = 0; My = 0;
    elseif step == 2
        Vx = 0; Vz = 1; My = 0;
    elseif step == 3
        Vx = 0; Vz = 0; My = 1;
    end
    
    %Calculate inertia products
    term2_P = (Ix * Vx - Ixz * Vz)/(Ix*Iz - Ixz^2);
    
    term4_P = (Iz * Vz - Ixz * Vx)/(Ix*Iz - Ixz^2);
    
    %term2_S = (Ix *
    
    %term4_S =
    
    %Section 1
    for i = 1:(num_str_1 + 2)
        term1_P = Section1.Areas(i) .* (Section1.Loc_Str(i,1) - XC); %A(x-xc)
        term3_P = Section1.Areas(i) .* (Section1.Loc_Str(i,2) - ZC); %A(z-zc)
        Section1.dP(i) = -term1_P * (term2_P) - term3_P * term4_P;
    end
    
    %Section 2
    for i = 1:(num_str_2 + 2)
        term1_P = Section2.Areas(i) .* (Section2.Loc_Str(i,1) - XC); %A(x-xc)
        term3_P = Section2.Areas(i) .* (Section2.Loc_Str(i,2) - ZC); %A(z-zc)
        Section2.dP(i) = -term1_P * (term2_P) - term3_P * term4_P;
    end
    
    %Section 3
    for i = 1:(num_str_3)
        term1_P = Section3.Areas(i) .* (Section3.Loc_Str(i,1) - XC); %A(x-xc)
        term3_P = Section3.Areas(i) .* (Section3.Loc_Str(i,2) - ZC); %A(z-zc)
        Section3.dP(i) = -term1_P * term2_P - term3_P * term4_P;
    end
    
    %Section 4
    for i = 1:(num_str_4)
        term1_P = Section4.Areas(i) .* (Section4.Loc_Str(i,1) - XC); %A(x-xc)
        term3_P = Section4.Areas(i) .* (Section4.Loc_Str(i,2) - ZC); %A(z-zc)
        Section4.dP(i) = -term1_P * term2_P - term3_P * term4_P;
    end
    
    
    % Analysis - Divide sections into webs
    
    x1 = sc_1; %For calculating area btwn 3 points, first pt on camber
    z1 = camber(index_sc_1);
    
    %Start with Section 1
    for i = 1:(num_str_1 + 1)
        Webs1.start(i,1) = Section1.Loc_Str(i,1);
        Webs1.start(i,2) = Section1.Loc_Str(i,2);
        Webs1.end(i,1) = Section1.Loc_Str(i+1,1);
        Webs1.end(i,2) = Section1.Loc_Str(i+1,2);
        
        ind_1 = find(abs(x - Webs1.start(i,1)) < 0.001);
        ind_2 = find(abs(x - Webs1.end(i,1)) < 0.001);
        vec_x = x(ind_1:ind_2);
        vec_z = z_upper(ind_1:ind_2);
        Webs1.len(i) = get_len(vec_x, vec_z);
        Webs1.t = Section1.Thickness;
        
        x2 = Webs1.start(i,1);
        z2 = Webs1.start(i,2);
        x3 = Webs1.end(i,1);
        z3 = Webs1.end(i,2);
        Webs1.Areas(i) = polyarea([x1 x2 x3], [z1 z2 z3]); %Swept area btwn 3 pts
        
        %For incremental area trapz(z1, z2-z1) where z1 is begin/end and z is
        %distance between the two curve
    end
    
    % scatter(x1, z1, 'bx');
    % scatter(x2, z2, 'bx'); %Is this point correct??
    % scatter(x3, z3, 'bx');
    
    %Section 2
    for i = 1:(num_str_2 + 1)
        Webs2.start(i,1) = Section2.Loc_Str(i,1);
        Webs2.start(i,2) = Section2.Loc_Str(i,2);
        Webs2.end(i,1) = Section2.Loc_Str(i+1,1);
        Webs2.end(i,2) = Section2.Loc_Str(i+1,2);
        
        ind_1 = find(abs(x - Webs2.start(i,1)) < 0.001);
        ind_2 = find(abs(x - Webs2.end(i,1)) < 0.001);
        vec_x = x(ind_1:ind_2);
        vec_z = z_lower(ind_1:ind_2);
        Webs2.len(i) = get_len(vec_x, vec_z);
        Webs2.t = Section2.Thickness;
        x2 = Webs2.start(i,1);
        z2 = Webs2.start(i,2);
        x3 = Webs2.end(i,1);
        z3 = Webs2.end(i,2);
        Webs2.Areas(i) = polyarea([x1 x2 x3], [z1 z2 z3]); %Swept area btwn 3 pts
    end
    
    %Section 3
    for i = 1:(num_str_3)
        Webs3.start(i,1) = Section3.Loc_Str(i,1);
        Webs3.start(i,2) = Section3.Loc_Str(i,2);
        if i == num_str_3
            Webs3.end(i,1) = Section2.Loc_Str(1,1);
            Webs3.end(i,2) = Section2.Loc_Str(1,2);
        else
            Webs3.end(i,1) = Section3.Loc_Str(i+1,1);
            Webs3.end(i,2) = Section3.Loc_Str(i+1,2);
        end
        ind_1 = find(abs(x - Webs3.start(i,1)) < 0.001);
        ind_2 = find(abs(x - Webs3.end(i,1)) < 0.001);
        vec_x = x(ind_1:ind_2);
        vec_z = z_lower(ind_1:ind_2);
        Webs3.len(i) = get_len(vec_x, vec_z);
        Webs3.t = Section3.Thickness;
        x2 = Webs3.start(i,1);
        z2 = Webs3.start(i,2);
        x3 = Webs3.end(i,1);
        z3 = Webs3.end(i,2);
        Webs3.Areas(i) = polyarea([x1 x2 x3], [z1 z2 z3]); %Swept area btwn 3 pts
    end
    
    %Section 4
    for i = 1:(num_str_4)
        Webs4.start(i,1) = Section4.Loc_Str(i,1);
        Webs4.start(i,2) = Section4.Loc_Str(i,2);
        if i == num_str_4
            Webs4.end(i,1) = Section1.Loc_Str(1,1);
            Webs4.end(i,2) = Section1.Loc_Str(1,2);
        else
            Webs4.end(i,1) = Section4.Loc_Str(i+1,1);
            Webs4.end(i,2) = Section4.Loc_Str(i+1,2);
        end
        ind_1 = find(abs(x - Webs4.start(i,1)) < 0.001);
        ind_2 = find(abs(x - Webs4.end(i,1)) < 0.001);
        vec_x = x(ind_1:ind_2);
        vec_z = z_upper(ind_1:ind_2);
        Webs4.len(i) = get_len(vec_x, vec_z);
        Webs4.t = Section4.Thickness;
        x2 = Webs4.start(i,1);
        z2 = Webs4.start(i,2);
        x3 = Webs4.end(i,1);
        z3 = Webs4.end(i,2);
        Webs4.Areas(i) = polyarea([x1 x2 x3], [z1 z2 z3]); %Swept area btwn 3 pts
    end
    
    %Left end
    WebLeft.start(1,1) = Section3.Loc_Str(1,1);
    WebLeft.start(1,2) = Section3.Loc_Str(1,2);
    WebLeft.end(1,1) = Section4.Loc_Str(1,1);
    WebLeft.end(1,2) = Section4.Loc_Str(1,2);
    
    ind_2_top = find(abs(x - WebLeft.end(1,1)) < 0.001);
    ind_1_bot = find(abs(x - WebLeft.start(1,1)) < 0.001);
    
    vec_x_bot = x(1:ind_1_bot);
    vec_z_bot = z_lower(1:ind_1_bot);
    vec_x_top = x(1:ind_2_top);
    vec_z_top = z_upper(1:ind_2_top);
    WebLeft.len = arclength(vec_x_bot, vec_z_bot) + arclength(vec_x_top, vec_z_top);
    WebLeft.Area = polyarea([x1 WebLeft.start(1,1) WebLeft.end(1,1)], [z1 WebLeft.start(1,2) WebLeft.end(1,2)])...
        + polyarea([0 WebLeft.start(1,1) WebLeft.end(1,1)], [0 WebLeft.start(1,2) WebLeft.end(1,2)]); %Two polyareas
    WebLeft.t = Section3.Thickness;
    
    %Right end
    WebRight.start(1,1) = Section1.Loc_Str(num_str_1 + 2,1);
    WebRight.start(1,2) = Section1.Loc_Str(num_str_1 + 2,2);
    WebRight.end(1,1) = Section2.Loc_Str(num_str_2 + 2,1);
    WebRight.end(1,2) = Section2.Loc_Str(num_str_2 + 2,2);
    WebRight.len = norm(WebRight.start - WebRight.end);
    WebRight.Area = polyarea([x1 WebRight.start(1,1) WebRight.end(1,1)], [z1 WebRight.start(1,2) WebRight.end(1,2)]);
    WebRight.t = SparCap.Thickness2;
    
    % Find shears
    %Got to do it for the left end and right end and then compile into a single
    %structure that runs clockwise (so Section 1, right end, reverse of Section
    %2, Reverse of Section 3, left end then Section 4)
    
    Per = [Webs1.start; WebRight.start; flipud(Webs2.end); flipud(Webs3.end); WebLeft.start; Webs4.start];
    Lens = transpose([Webs1.len WebRight.len fliplr(Webs2.len) fliplr(Webs3.len) WebLeft.len Webs4.len]);
    Area = transpose([Webs1.Areas WebRight.Area fliplr(Webs2.Areas) fliplr(Webs3.Areas) WebLeft.Area Webs4.Areas]);
    
    dP_right = [transpose(Section1.dP); flipud(transpose(Section2.dP))];
    dP_left = vertcat(Section2.dP(1), flipud(transpose(Section3.dP)), transpose(Section4.dP), Section1.dP(1));
    
    Web = struct('Length', Lens, 'Area', Area, 'Coords', Per, 'P_right', dP_right, 'P_left', dP_left);
    
    %Calculate q' values
    
    %Right
    Web.q_right(1) = 0; %We will cut the first web starting from spar cap 1
    
    for i = 2:(num_str_1 + num_str_2 + 4)
        Web.q_right(i) = Web.q_right(i-1) + Web.P_right(i);
    end
    
    Web.q_left(1) = 0;
    
    for i = 2:(num_str_3 + num_str_4 + 2)
        Web.q_left(i) = Web.q_left(i-1) + Web.P_left(i);
    end
    
    Web.q_right = transpose(Web.q_right);
    Web.q_left = transpose(Web.q_left);
    % Calculate q * ds/t
    
  
    %Section 1
    for i = 1:(num_str_1 + 1)
        ds_t_1(i) = Web.Length(i) ./  Section1.Thickness;
    end
    
    %Right
    ds_t_R = Web.Length(num_str_1 + 2) ./ SparCap.Thickness2;
    
    %Section 2
    for j = 1:(num_str_2 + 1)
        ds_t_2(j) = Web.Length(num_str_1 + 2 + j) ./ Section2.Thickness;
    end
    
    % %Mid Web
    WebMid.start(1,1) = Section1.Loc_Str(1,1);
    WebMid.start(1,2) = Section1.Loc_Str(1,2);
    WebMid.end(1,1) = Section2.Loc_Str(1,1);
    WebMid.end(1,2) = Section2.Loc_Str(1,2);
    WebMid.len = norm(WebMid.start - WebMid.end);
    ds_t_M = WebMid.len / SparCap.Thickness1;
    
    %Section 3
    for k = 1:(num_str_3)
        ds_t_3(k) = Web.Length(num_str_1 + 2 + num_str_2 + 1 + k) ./ Section3.Thickness;
    end
    
    %Left
    ds_t_L = Web.Length(num_str_1 + 2 + num_str_2 + 1 + num_str_3 + 1)./Section3.Thickness;
    %Section 4
    for l = 1:(num_str_4)
        ds_t_4(l) = Web.Length(num_str_1 + 2 + num_str_2 + 1 + num_str_3 + 1 + l) ./ Section4.Thickness;
    end
    
    Web.ds_t = vertcat(transpose(ds_t_1), ds_t_R, transpose(ds_t_2), transpose(ds_t_3), ds_t_L, transpose(ds_t_4), ds_t_M);
    %ds_t end is the spar 1
    
    % Calculate q * ds/t

    %For the right
    for i = 1:length(Web.P_right)
        if i == length(Web.P_right)
            qdst_r(i) = Web.q_right(i) .* Web.ds_t(end);
        else
            qdst_r(i) =  Web.q_right(i) .* Web.ds_t(i);
        end
    end
    
    Web.qdst_r = transpose(qdst_r);
    
    %For the left
    for i = 1:length(Web.P_left)
        if i == length(Web.P_left)
            qdst_l(i) = Web.q_left(i) .* Web.ds_t(end);
        else
            qdst_l(i) =  Web.q_left(i) .* Web.ds_t(i+length(Web.P_right)-1);
        end
    end
    
    Web.qdst_l = transpose(qdst_l);
    
    %Summing up all the q ds/t
    Web.rqsum = sum(Web.qdst_r);
    Web.lqsum = sum(Web.qdst_l);
    
    %Summing up all the ds/t
    Web.rsum = sum(Web.ds_t(1:length(Web.P_right)-1)) + Web.ds_t(end);
    Web.lsum = sum(Web.ds_t(length(Web.P_right):end));
    
    % solve this mother fucker
    
    A = [Web.rsum -Web.ds_t(end);-Web.ds_t(end) Web.lsum];
    B = [-Web.rqsum; -Web.lqsum];
    qs = (A\B);
    
    % Moment balance
    
    %First need to find area of left and right cell, starting with right cell
    Z_left_top = z_upper(1:index_sc_1);
    Z_left_lower = z_lower(1:index_sc_1);
    X_left = x(1:index_sc_1);
    Area_left = trapz(X_left, Z_left_top) + abs(trapz(X_left, Z_left_lower));
    
    %Check area left
    Area_left2 = sum(Web.Area(length(Web.q_right):end));
    
    Z_right_top = z_upper(index_sc_1:index_sc_2);
    Z_right_lower = z_lower(index_sc_1:index_sc_2);
    X_right = x(index_sc_1:index_sc_2);
    Area_right = trapz(X_right, Z_right_top) + abs(trapz(X_right, Z_right_lower));
    
    %Check area right
    Area_right2 = sum(Web.Area(1:length(Web.q_right)));

    %Calculate the torques
    torque_qs = 2 * Area_right * qs(1) + 2 * Area_left * qs(2);
   
    torque_qp_right = (2 *(Web.q_right) .* Web.Area(1:(num_str_1 + num_str_2 + 4)));
    torque_qp_right(end) = 0;
    torque_qp_left = (2 .*(Web.q_left(1:num_str_3+num_str_4+1)) .* Web.Area((num_str_1 + num_str_2 + 4:end)));
%     torque_qp_left(end) = 0;
    
    torque_qp = sum(torque_qp_right) + sum(torque_qp_left);
    torque = torque_qs + torque_qp;

    if step == 1
        shc(step) = torque/Vx;
        shc_1(2) = z1 + shc(step); %normalize against location of area sweeps (x1, z1)%z location
    else
        shc(step) = torque/Vz;
        shc_1(1) = x1 - shc(step);  %shc(2) is x location of shear center
    end
    if step == 2
        shearplot = scatter(shc_1(1), shc_1(2), 'kx');
        xloc = shearplot.XData; %For future use
    end
  
    
    %Find q_t
    
    if step == 1
        torsion(step) = Vx * (shc_1(step));
    elseif step == 2
        torsion(step) = Vz * (shc_1(step) - 1.25);
    elseif step == 3
        torsion(step) = My;
    end
    
    %Solving the twist equation on the right/left cells
    t_right = Web.rsum / Area_right;
    t_left = Web.lsum / Area_left;
    t_left_mid = Web.ds_t(end) / Area_left;
    t_right_mid = Web.ds_t(end) / Area_right;
    q1t_q2t = (t_right + t_right_mid) / (t_left + t_left_mid);
    

    q2t = torsion(step) / (2 * Area_left * q1t_q2t + 2 * Area_right);
    q1t = q2t * q1t_q2t;
    
    %Getting total shear flow
    for i = 1:length(Web.q_right)
        Web.QR(i) = Web.q_right(i) + q1t + qs(1);
    end
    
    for i = 1:length(Web.q_left)
        Web.QL(i) = Web.q_left(i) + q2t + qs(2);
    end
    
    Web.QR = transpose(Web.QR);
    Web.QL = transpose(Web.QL);
    q_unit(:,step) = [Web.QR(1:end-1); Web.QL];
end
  
hold off;

%% Calculate bending stresses and column buckling

E = 1439006421.41; %lbs/ft^2, Young's modulus converted from 69 GPa
k = 1.5; %Based on supports

numRibs = randi(15); %number of ribs = random integer from 1:8
numRibs = 11;
distRibs = 17/(numRibs + 1); %Distance between ribs, ft

locRibs(1) = 0;
for i = 1:(numRibs + 1)
    locRibs(i+1) = locRibs(i) + distRibs;
end

%Find the index at which each rib is at to find the corresponding Vs, Ms
y = round(Load.L1.y, 1);
z = round(locRibs, 1);
for i = 1:(numRibs + 2)
    a = y - round(locRibs(i),1);
    indexRibs(i) = find(a==min(abs(a)),1);
end

%Import loads

for i = 1:length(indexRibs)
    Vx_3(i) = Load.L3.Vx(indexRibs(i)); 
    Vz_3(i) = Load.L3.Vz(indexRibs(i));
    Mz_3(i) = Load.L3.Mz(indexRibs(i));
    My_3(i) = Load.L3.My(indexRibs(i));
    Mx_3(i) = Load.L3.Mx(indexRibs(i)); %This has the largest!!
end

%Calculate out your bending stresses
% term2_S_1 = (Ix * Mz_1 + Ixz * Mx_1)/(Ix*Iz - Ixz^2);
% term2_S_2 = (Ix * Mz_2 + Ixz * Mx_2)/(Ix*Iz - Ixz^2);
term2_S_3 = (Ix * Mz_3 + Ixz * Mx_3)/(Ix*Iz - Ixz^2);
% term2_S_4 = (Ix * Mz_4 + Ixz * Mx_4)/(Ix*Iz - Ixz^2);
% term2_S_5 = (Ix * Mz_5 + Ixz * Mx_5)/(Ix*Iz - Ixz^2);
% term2_S_6 = (Ix * Mz_6 + Ixz * Mx_6)/(Ix*Iz - Ixz^2);

% term4_S_1 = (Iz * Mx_1 + Ixz * Mz_1)/(Ix*Iz - Ixz^2);
% term4_S_2 = (Iz * Mx_2 + Ixz * Mz_2)/(Ix*Iz - Ixz^2);
term4_S_3 = (Iz * Mx_3 + Ixz * Mz_3)/(Ix*Iz - Ixz^2);
% term4_S_4 = (Iz * Mx_4 + Ixz * Mz_4)/(Ix*Iz - Ixz^2);
% term4_S_5 = (Iz * Mx_5 + Ixz * Mz_5)/(Ix*Iz - Ixz^2);
% term4_S_6 = (Iz * Mx_6 + Ixz * Mz_6)/(Ix*Iz - Ixz^2);

%Note that Section.sigy contains a numRibsx6xnumStr containing the bending
%stresses

%Section 1
for i = 1:(num_str_1 + 2)
    term1_S = Section1.Loc_Str(i,1) - XC; %(x-xc)
    term3_S = Section1.Loc_Str(i,2) - ZC; %(z-zc)
    Section1.sigy_3(:,i) = transpose(term1_S * (term2_S_3) - term3_S * term4_S_3);
    Section1.Ix(i) = Section1.Areas(i) * (term3_S)^2;
    Section1.Iz(i) = Section1.Areas(i) * (term1_S)^2;
    Section1.Fcr(i) = (pi^2 * E * Section1.Ix(i))/((k*distRibs)^2);
    Section1.sigccr(i) = Section1.Fcr(i) ./ Section1.Areas(i);
end

%Section 2
for i = 1:(num_str_2 + 2)
    term1_S = Section2.Loc_Str(i,1) - XC; %(x-xc)
    term3_S = Section2.Loc_Str(i,2) - ZC; %(z-zc)
    Section2.sigy_3(:,i) = transpose(term1_S * (term2_S_3) - term3_S * term4_S_3);
    Section2.Ix(i) = Section2.Areas(i) * (term3_S)^2;
    Section2.Iz(i) = Section2.Areas(i) * (term1_S)^2;
    Section2.Fcr(i) = (pi^2 * E * Section2.Ix(i))/((k*distRibs)^2);
    Section2.sigccr(i) = Section2.Fcr(i) ./ Section2.Areas(i);
end

%Section 3
for i = 1:(num_str_3)
    term1_S = Section3.Loc_Str(i,1) - XC; %(x-xc)
    term3_S = Section3.Loc_Str(i,2) - ZC; %(z-zc)
    Section3.sigy_3(:,i) = transpose(term1_S * (term2_S_3) - term3_S * term4_S_3);
    Section3.Ix(i) = Section3.Areas(i) * (term3_S)^2;
    Section3.Iz(i) = Section3.Areas(i) * (term1_S)^2;
    Section3.Fcr(i) = (pi^2 * E * Section3.Ix(i))/((k*distRibs)^2);
    Section3.sigccr(i) = Section3.Fcr(i) ./ Section3.Areas(i);
end

%Section 4
for i = 1:(num_str_4)
    term1_S = Section4.Loc_Str(i,1) - XC; %(x-xc)
    term3_S = Section4.Loc_Str(i,2) - ZC; %(z-zc)
    Section4.sigy_3(:,i) = transpose(term1_S * (term2_S_3) - term3_S * term4_S_3);
    Section4.Ix(i) = Section4.Areas(i) * (term3_S)^2;
    Section4.Iz(i) = Section4.Areas(i) * (term1_S)^2;
    Section4.Fcr(i) = (pi^2 * E * Section4.Ix(i))/((k*distRibs)^2);
    Section4.sigccr(i) = Section4.Fcr(i) ./ Section4.Areas(i);
end

%Compare against sig_y

%Calculate out your average sigy for each rib web and FOS
for j = 1:size(Section1.sigy_3,2)
    for i = 1:(numRibs + 1)
        Section1.sigave(i,j) = (Section1.sigy_3(i,j) + Section1.sigy_3(i+1,j)) / 2;
        Section1.FOS(i,j) = Section1.sigave(i,j) ./ Section1.sigccr(j);
        Section1.FOS(i,j) = 1/Section1.FOS(i,j);
    end
end

for j = 1:size(Section2.sigy_3,2)
    for i = 1:(numRibs + 1)
        Section2.sigave(i,j) = (Section2.sigy_3(i,j) + Section2.sigy_3(i+1,j)) / 2;
        Section2.FOS(i,j) = Section2.sigave(i,j) ./ Section2.sigccr(j);
        Section2.FOS(i,j) = 1/Section2.FOS(i,j);
    end
end

for j = 1:size(Section3.sigy_3,2)
    for i = 1:(numRibs + 1)
        Section3.sigave(i,j) = (Section3.sigy_3(i,j) + Section3.sigy_3(i+1,j)) / 2;
        Section3.FOS(i,j) = Section3.sigave(i,j) ./ Section3.sigccr(j);
        Section3.FOS(i,j) = 1/Section3.FOS(i,j);
    end
end

for j = 1:size(Section4.sigy_3,2)
    for i = 1:(numRibs + 1)
        Section4.sigave(i,j) = (Section4.sigy_3(i,j) + Section4.sigy_3(i+1,j)) / 2;
        Section4.FOS(i,j) = Section4.sigave(i,j) ./ Section4.sigccr(j);
        Section4.FOS(i,j) = 1/Section4.FOS(i,j);
    end
end

%% Panel buckling
%Need to apply across ribs
a_b = 2; %For future reference
k_comp = 7.75;
k_bend = 36;
k_shear = 8;
v = 0.33;
shear = q_unit(:,1) * Vx_3 + q_unit(:,2) * Vz_3 + q_unit(:,3) * My_3;
shear = transpose(shear); % [numRibs x numStringers] matrix with total shear values
%Calculate out your total shear forces

%Section 1
for j = 1:(numRibs + 2)
    for i = 1:(Section1.Num_Str - 1)
        Section1.panel_shear(j,i) = shear(j,i) / Section1.Thickness;
    end
end

for j = 1:(numRibs + 1)
    for i = 1:(Section1.Num_Str - 1)
        Section1.panel_ave(j,i) = (Section1.panel_shear(j,i) + Section1.panel_shear(j+1,i)) / 2; %Each panel in Section 1 (including ribs) average stress
        
        Section1.aveComp(j,i) = (Section1.sigave(j,i) + Section1.sigave(j,i+1))/2; %Obtain 4 compressive stress from based on 4 stresses
        Section1.sigpcr_comp(j,i) = (k_comp * pi^2 * E)/(12*(1 - v^2)) * (Section1.Thickness/Web.Length(i))^2;
        Section1.sigpcr_bend(j,i) = (k_bend * pi^2 * E)/(12*(1 - v^2)) * (Section1.Thickness/Web.Length(i))^2;
        Section1.sigpcr_shear(j,i) = (k_shear * pi^2 * E)/(12*(1 - v^2)) * (Section1.Thickness/Web.Length(i))^2;
        
        Section1.FOS_comp(j,i) = Section1.sigpcr_comp(j,i) ./ Section1.aveComp(j,i);
        Section1.FOS_shear(j,i) = Section1.sigpcr_shear(j,i) ./ Section1.panel_ave(j,i);
    end
end

%Section 2
incr = num_str_1 + 2; %Used to adjust index for shear
for j = 1:(numRibs + 2)
    for i = 1:(Section2.Num_Str - 1)
        Section2.panel_shear(j,i) = shear(j,i + incr) / Section2.Thickness;
    end
end

for j = 1:(numRibs + 1)
    for i = 1:(Section2.Num_Str - 1)
        Section2.panel_ave(j,i) = (Section2.panel_shear(j,i) + Section2.panel_shear(j+1,i)) / 2; %Each panel in Section 1 (including ribs) average stress
        
        Section2.aveComp(j,i) = (Section2.sigave(j,i) + Section2.sigave(j,i+1))/2; %Obtain 4 compressive stress from based on 4 stresses
        Section2.sigpcr_comp(j,i) = (k_comp * pi^2 * E)/(12*(1 - v^2)) * (Section2.Thickness/Web.Length(i + incr))^2;
        Section2.sigpcr_bend(j,i) = (k_bend * pi^2 * E)/(12*(1 - v^2)) * (Section2.Thickness/Web.Length(i + incr))^2;
        Section2.sigpcr_shear(j,i) = (k_shear * pi^2 * E)/(12*(1 - v^2)) * (Section2.Thickness/Web.Length(i + incr))^2;
        
        Section2.FOS_comp(j,i) = Section2.sigpcr_comp(j,i) ./ Section2.aveComp(j,i);
        Section2.FOS_shear(j,i) = Section2.sigpcr_shear(j,i) ./ Section2.panel_ave(j,i);
    end
end

%Section 3
incr = num_str_1 + num_str_2 + 3; %Used to adjust index for shear
for j = 1:(numRibs + 2)
    for i = 1:(Section3.Num_Str) %Check this
        Section3.panel_shear(j,i) = shear(j,i + incr) / Section3.Thickness;
    end
end

for j = 1:(numRibs + 1)
    for i = 1:(Section3.Num_Str - 1)
        Section3.panel_ave(j,i) = (Section3.panel_shear(j,i) + Section3.panel_shear(j+1,i)) / 2; %Each panel in Section 1 (including ribs) average stress
        
        Section3.aveComp(j,i) = (Section3.sigave(j,i) + Section3.sigave(j,i+1))/2; %Obtain 4 compressive stress from based on 4 stresses
        Section3.sigpcr_comp(j,i) = (k_comp * pi^2 * E)/(12*(1 - v^2)) * (Section3.Thickness/Web.Length(i + incr))^2;
        Section3.sigpcr_bend(j,i) = (k_bend * pi^2 * E)/(12*(1 - v^2)) * (Section3.Thickness/Web.Length(i + incr))^2;
        Section3.sigpcr_shear(j,i) = (k_shear * pi^2 * E)/(12*(1 - v^2)) * (Section3.Thickness/Web.Length(i + incr))^2;
        
        Section3.FOS_comp(j,i) = Section3.sigpcr_comp(j,i) ./ Section3.aveComp(j,i);
        Section3.FOS_shear(j,i) = Section3.sigpcr_shear(j,i) ./ Section3.panel_ave(j,i);
    end
end

%Section 4
incr = num_str_1 + num_str_2 + num_str_3 + 3; %Used to adjust index for shear
for j = 1:(numRibs + 2)
    for i = 1:(Section4.Num_Str - 1) %Check this
        Section4.panel_shear(j,i) = shear(j,i + incr) / Section4.Thickness;
    end
end

for j = 1:(numRibs + 1)
    for i = 1:(Section4.Num_Str - 1)
        Section4.panel_ave(j,i) = (Section4.panel_shear(j,i) + Section4.panel_shear(j+1,i)) / 2; %Each panel in Section 1 (including ribs) average stress
        
        Section4.aveComp(j,i) = (Section4.sigave(j,i) + Section4.sigave(j,i+1))/2; %Obtain 4 compressive stress from based on 4 stresses
        Section4.sigpcr_comp(j,i) = (k_comp * pi^2 * E)/(12*(1 - v^2)) * (Section4.Thickness/Web.Length(i + incr))^2;
        Section4.sigpcr_bend(j,i) = (k_bend * pi^2 * E)/(12*(1 - v^2)) * (Section4.Thickness/Web.Length(i + incr))^2;
        Section4.sigpcr_shear(j,i) = (k_shear * pi^2 * E)/(12*(1 - v^2)) * (Section4.Thickness/Web.Length(i + incr))^2;
        
        Section4.FOS_comp(j,i) = Section4.sigpcr_comp(j,i) ./ Section4.aveComp(j,i);
        Section4.FOS_shear(j,i) = Section4.sigpcr_shear(j,i) ./ Section4.panel_ave(j,i);
    end
end

%% Heatmap the FOS's and find minimum FOS
minFOS = 0;

figure
heatmap(abs([Section4.FOS_comp Section1.FOS_comp]), 'Along Cross Section','Span', 1);
title('Distribution of failure margins for compressive failure along top of wing');

figure
heatmap(abs([Section3.FOS_comp Section2.FOS_comp]), 'Along Cross Section','Span', 1);
title('Distribution of failure margins for compressive failure along bottom of wing');

minFOS = min(min(abs(Section1.FOS_comp)));

figure
heatmap(abs([Section4.FOS_shear Section1.FOS_shear]), 'Along Cross Section','Span', 1);
title('Distribution of failure margins for shear failure along top of wing');

figure
heatmap(abs([Section3.FOS_shear Section2.FOS_shear]), 'Along Cross Section','Span', 1);
title('Distribution of failure margins for shear failure along bottom of wing');

temp = min(min(abs(Section1.FOS_shear)));
if temp < minFOS
    minFOS = temp;
end

figure
heatmap(abs([Section4.FOS Section1.FOS]), 'Along Cross Section','Span', 1);
title('Distribution of failure margins for bending failure along top of wing');

figure
heatmap(abs([Section3.FOS Section2.FOS]), 'Along Cross Section','Span', 1);
title('Distribution of failure margins for bending failure along bottom of wing');

temp = min(min(abs(Section1.FOS)));
if temp < minFOS
    minFOS = temp;
end


%% Wing Divergence
b_2 = 17;
G = 543021391; %lb/ft^2

%Right
A2R = 4 * Area_right^2;
dstr = sum(Web.ds_t(1:num_str_1 + num_str_2 + 3)) + Web.ds_t(end);
GJR = A2R * G / dstr;

%Left
A2L = 4 * Area_left^2;
dstl = sum(Web.ds_t(num_str_1 + num_str_2 + 4:end));
GJL = A2L * G / dstl;

%Divergence
GJ = GJR + GJL;
rho_d = 0.002238; %slugs/ft^3
CL_a = 6.73; %From master code
e = abs(xloc - 0.25 * c)/c;
num = pi^2 * GJ;
den = 2 * rho_d * CL_a * e * c^2 * (b_2)^2;
V_diverge = sqrt(num/den);

%% Weight
%Mass of the outer perimeter
total_length = sum([Web.Length]);
average_t = (Section1.Thickness + Section2.Thickness + Section3.Thickness + Section4.Thickness)/4;
Volume = 17 * total_length * average_t; %ft^3
Density = 5.24; %slugs/ft^3
Mass = Volume * Density; %slug

%Mass of stringers
Volume_str = (sum(Section1.Areas) + sum(Section2.Areas) + sum(Section3.Areas) + sum(Section4.Areas)) * 17;
Mass = Mass + Density * Volume_str;

%Mass of ribs
Mass = Mass + 15 * numRibs;

%% FINAL SPEC OUTPUT FOR MONTE CARLO ANALYSIS

Wing.Section1 = Section1;
Wing.Section2 = Section2;
Wing.Section3 = Section3;
Wing.Section4 = Section4;
Wing.VD = V_diverge; 
Wing.Weight = Mass;
Wing.FOS = minFOS;

%% MONTE CARLO OPTIMIZATION ROUTINE
% %% MONTE CARLO ROUTINE
% iterations = 1;
% tic
% 
% for mc = 1:iterations
%     Structural_Analysis_Master_v6
%     it = ['It_',num2str(mc)];
%     Monte.(it) = Wing; %Saves all results into this structure
%     W(mc) = Wing.Weight;
%     FOS(mc) = Wing.FOS;
% end
% end_time = toc;
% X = ['Number of iterations: ', num2str(iterations)];
% disp(X);
% X = ['Average time per iteration: ',num2str(end_time/iterations), ' seconds'];
% disp(X);
% 
% % %Plot Weights
% % % subplot(1,2,1);
% % hist(W);
% % X = ['Distribution of weights for ', num2str(iterations), ' simulations'];
% % title(X);
% % xlabel('Weight of wing (lbs)');
% % ylabel('Frequency');
% % 
% % %Scatterplot your FOS against Weight
% % % subplot(1,2,2);
% % figure
% % scatter(W, FOS, 'filled');
% % xlabel('Weights (lbs)');
% % ylabel('Minimum margins');
% % X = ['Scatterplot of minimum margins of failure against weights for ', num2str(iterations), ' simulations'];
% % title(X);
