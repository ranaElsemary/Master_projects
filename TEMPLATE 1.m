%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following script executes a VLM to study the aerodynamic properties
% of a finite swept wing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all

%% Figure's positions
% ----------------------------------------------------------------------- %
scrsz = get(0,'ScreenSize');
scW=scrsz(3); % Screen Width (px)
scH=scrsz(4); % Screen Height (px)
figure('Name','Wing','Position',[1 scH/4 scW/2 scH/2.5])
figure('Name','Distributions','Position',[scW/2 scH/4 scW/2 scH/2.5])

% Plots configuration
font=15;  % font size
lw=2;     % linewidth
szp=100;  % scatter size plot

%% Mesh Configuration
% ----------------------------------------------------------------------- %
npx = 2; % Number of panels in the streamwise direction
npy = 8; % Number of panels in the SEMI-spanwise direction
Num_panels = npx*npy;

%% Flight Conditions
% ----------------------------------------------------------------------- %

alpha = 1*pi/180; % Angle of Attack (rad)
Uinf  = 20;       % Wind Speed (m/s)
rho   = 1.225;    % Air Density (kg/m3)

%% Wing Geometry
% ----------------------------------------------------------------------- %

Ale = deg2rad(20);            % Leading edge sweep angle (rad)
phi = deg2rad(10)+ (1e-10);  % Dihedral angle (rad)
s   = 11.43/2;            % SEMI-span (m)
cr  = 4.04;            % Root chord (m)
ct  = 1.68 + (1e-10);  % Tip chord (m)
S   = 37.16;            % Wing surface (m)
AR  = 3.52;            % Wing Aspect Ratio (-)

% IMPORTANT: (1e-10) added at Dihedral Angle and Tip Chord in order to
% avoid numerical issues.

%% Panels´ Geometrical Properties
% ----------------------------------------------------------------------- %

tic %Initializes time measuring
for j = 1:npy
    
    % Local y coordinate of each slice 
    y_A_xy= (s/npy)*(j-1);
    y_B_xy= (s/npy)*(j);
    % Local chord at the coordinate y 
    localChord_A=cr-(y_A_xy*cr/s)-(y_A_xy*ct/s);
    localChord_B=cr-(y_B_xy*cr/s)-(y_B_xy*ct/s);

    % Local dx
    dx_A(j)=localChord_A/npx; 
    dx_B(j)=localChord_B/npx; 

    dx_A(j) = cr/npx-(j-1)*(cr-ct)/(npx*npy);
    dx_B(j) = cr/npx-j*(cr-ct)/(npx*npy);

    for i = 1:npx      
    % Point A of each panel [y1n,z1n] -> End A of each bound vortex filament
    Panels(i,j).y1n = y_A_xy;
    Panels(i,j).z1n = abs(y_A_xy)*tan(phi);   

    % Point B of each panel [y2n,z2n]-> End B of each vortex filament
    Panels(i,j).y2n = Panels(i,j).y1n + (s/npy);
    Panels(i,j).z2n = Panels(i,j).z1n + (s/npy)*tan(phi);     
    
    % Control Points Locations [ym,zm]
    Panels(i,j).ym = (Panels(i,j).y1n + Panels(i,j).y2n)/2 ;
    Panels(i,j).zm = (Panels(i,j).z1n + Panels(i,j).z2n) /2 ; 
    end
    % Control Points Locations [x1n,x2n,xm]
    for i = 1:npx
        Panels(i,j).x1n = 0.25 *dx_A(j) + dx_A(j)*(i-1)+Panels(i,j).y1n*tan(Ale);
        Panels(i,j).x2n = 0.25 *dx_B(j) + dx_B(j)*(i-1)+Panels(i,j).y2n * tan(Ale);
        Panels(i,j).xm = 0.375 *(dx_A(j)+dx_B(j)) + 0.5*(dx_A(j)+dx_B(j))*(i-1)+Panels(i,j).ym * tan(Ale);
    end
    end
     
    


% Time enlapsed in the geometric loop
T_geom_loop=toc; % (s)  


step = 0;
for n = 1:npy
    for m=1:npx
        step = step +1;
        %counter for components of A points 
        A_x(step) = Panels(m,n).x1n;
        A_y(step) = Panels(m,n).y1n;
        A_z(step) = Panels(m,n).z1n;
        %counter for components of B points
        B_x(step) = Panels(m,n).x2n;
        B_y(step) = Panels(m,n).y2n;
        B_z(step) = Panels(m,n).z2n;
        %counter for components of controll points
        C_x(step) = Panels(m,n).xm;
        C_y(step) = Panels(m,n).ym;
        C_z(step) = Panels(m,n).zm;
    end
end

%ploting the wing 
figure(1)
scatter3 (A_x, A_y, A_z, 'y', 'filled')
hold on
scatter3 (B_x, B_y, B_z, 'y', 'filled')
scatter3 (C_x, C_y, C_z, 'blue', 'filled')
scatter3 (A_x, -A_y, A_z, 'y', 'filled')
scatter3 (B_x, -B_y, B_z, 'y', 'filled')
scatter3 (C_x, -C_y, C_z, 'blue', 'filled')
for j = 1:npy
    dx = [Panels(1,j).x2n-0.25*dx_B(j),Panels(1,j).y2n,Panels(1,j).z2n;Panels(npx,j).x2n+0.75*dx_B(j) ,Panels(1,j).y2n,Panels(1,j).z2n];
    d_x = [Panels(1,j).x2n-0.25*dx_B(j),-Panels(1,j).y2n,Panels(1,j).z2n;Panels(npx,j).x2n+0.75*dx_B(j) ,-Panels(1,j).y2n,Panels(1,j).z2n];
    plot3 (dx(:,1),dx(:,2),dx(:,3),'b')
    plot3 (d_x(:,1),d_x(:,2),d_x(:,3),'b')
end
for i = 1:npx
    dy = [Panels(i,1).x1n-0.25*dx_A(1),Panels(i,1).y1n,Panels(i,1).z1n;Panels(i,npy).x2n-0.25*dx_B(npy) ,Panels(i,npy).y2n,Panels(i,npy).z2n];
    d_y = [Panels(i,1).x1n-0.25*dx_A(1),-Panels(i,1).y1n,Panels(i,1).z1n;Panels(i,npy).x2n-0.25*dx_B(npy) ,-Panels(i,npy).y2n,Panels(i,npy).z2n];
    plot3 (dy(:,1),dy(:,2),dy(:,3),'b')
    plot3 (d_y(:,1),d_y(:,2),d_y(:,3),'b')
end
    dx = [Panels(1,1).x1n-0.25*dx_A(1),Panels(1,1).y1n,Panels(1,1).z1n;Panels(npx,1).x1n+0.75*dx_A(1) ,Panels(1,1).y1n,Panels(1,1).z1n];    
    dy = [Panels(npx,1).x1n+0.75*dx_A(1),Panels(npx,1).y1n,Panels(npx,1).z1n;Panels(npx,npy).x2n+0.75*dx_B(npy) ,Panels(npx,npy).y2n,Panels(npx,npy).z2n];
    d_y = [Panels(npx,1).x1n+0.75*dx_A(1),-Panels(npx,1).y1n,Panels(npx,1).z1n;Panels(npx,npy).x2n+0.75*dx_B(npy) ,-Panels(npx,npy).y2n,Panels(npx,npy).z2n];
    plot3 (dx(:,1),dx(:,2),dx(:,3),'b')
    plot3 (dy(:,1),dy(:,2),dy(:,3),'b')
    plot3 (d_y(:,1),d_y(:,2),d_y(:,3),'b')
    pbaspect ([10,30,4]);

hold off

%% MAIN LOOP

% Set waitbar
h = waitbar(0,'Initializing Solver...');
% Reset Panel Counter
Panel=0;
tic
for i = 1:npx   
    for j = 1:npy % Solve only the Starboard wing due to symmetry
        % Panel Counter
        Panel=Panel+1;
        
        % Progress
        waitbar(Panel/(npx*npy),h,sprintf('Calculating %4.2f %%...',Panel/(npx*npy)*100))
        
        % Control point where the downwash from all the vortex is calculated
     
        xm =  Panels(i,j).xm;
        ym = Panels(i,j).ym ;
        zm =Panels(i,j).zm; 

        %% Starboard/Right Wing Vortex
        
        % Reset Vortex Counter
        Vortex=0; 
        
       for ii =1:npx %n
            for jj = 1: npy %m %right wing
                
                % Vortex Counter
                Vortex=Vortex+1;
                %Induced Velocity from the Bound Vortex
        
   A_z= [Panels(ii,jj).x1n Panels(ii,jj).y1n Panels(ii,jj).z1n] ;    % Point A of each panel         
   B= [Panels(ii,jj).x2n Panels(ii,jj).y2n Panels(ii,jj).z2n] ;    % Point B of each panel  
   C=  [Panels(i,j).xm, Panels(i,j).ym, Panels(i,j).zm];        % control points
   D=  [100000000000, 0, 0]; % infinty point
   r0= B-A_z; %AB
   r1= C-A_z; %AC
   r2= C-B; %BC
   r3= A_z-D; %DA
   r4= C-D; %DC
   r5= D-B; %BD
   r6= D-C; %CD
   fac1_AB= cross(r1,r2)/norm(cross(r1,r2))^2;
   fac2_AB= (dot(r0,r1)/norm(r1) - dot(r0, r2)/norm(r2));
   fac1_AD= cross(r4,r1)/norm(cross(r4,r1))^2;
   fac2_AD= (dot(r3,r4)/norm(r4) - dot(r3, r1)/norm(r1));
   fac1_BD= cross(r2,r4)/norm(cross(r2,r4))^2;
   fac2_BD= (dot(r5,r2)/norm(r2) - dot(r5, r4)/norm(r4));
   
   %The velocity induced by the whole horse shoe vortex for S wing
   Vi_s{i,j}(Vortex,:)= ((fac1_AB*fac2_AB)+(fac1_AD*fac2_AD)+(fac1_BD*fac2_BD))/(4*pi);

            end
        end
        
        %% Port/Left Wing Vortex
        
        % Reset Vortex Counter
        Vortex=0; 
        
         for ii =1:npx %n
            for jj = 1: npy %m %lift wing


                % Vortex Counter
                Vortex=Vortex+1;
                
            
       %Induced Velocity from the Bound Vortex
        
   A_z= [Panels(ii,jj).x1n -Panels(ii,jj).y1n Panels(ii,jj).z1n] ;    % Point A of each panel         
   B= [Panels(ii,jj).x2n -Panels(ii,jj).y2n Panels(ii,jj).z2n] ;    % Point B of each panel  
   C=  [Panels(i,j).xm, -Panels(i,j).ym, Panels(i,j).zm];        % control points
   D=  [100000000000, 0, 0]; % infinty point
   r0= B-A_z; %AB
   r1= C-A_z; %AC
   r2= C-B; %BC
   r3= A_z-D; %DA
   r4= C-D; %DC
   r5= D-B; %BD
   r6= D-C; %CD
   fac1_AB= cross(r1,r2)/norm(cross(r1,r2))^2;
   fac2_AB= (dot(r0,r1)/norm(r1) - dot(r0, r2)/norm(r2));
   fac1_AD= cross(r4,r1)/norm(cross(r4,r1))^2;
   fac2_AD= (dot(r3,r4)/norm(r4) - dot(r3, r1)/norm(r1));
   fac1_BD= cross(r2,r4)/norm(cross(r2,r4))^2;
   fac2_BD= (dot(r5,r2)/norm(r2) - dot(r5, r4)/norm(r4));
   
   %The velocity induced by the whole horse shoe vortex for S wing
   Vi_p{i,j}(Vortex,:)= ((fac1_AB*fac2_AB)+(fac1_AD*fac2_AD)+(fac1_BD*fac2_BD))/(4*pi);

 
                
            end
        end
           
    end
         
end

% Time enlapsed in the MAIN loop
T_MAIN_loop=toc; % (s) 

% Close wait bar
close(h)

%% Combine contributions from the Port and Startboard wings
Vortex=0;
for i=1:npx
    for j=1:npy
        Vortex=Vortex+1;
U(Vortex,:) = Vi_s{i,j}(:,1)+ Vi_p{i,j}(:,1);
V(Vortex,:) = Vi_s{i,j}(:,2)+Vi_p{i,j}(:,2);
W(Vortex,:) = Vi_s{i,j}(:,3)+ Vi_p{i,j}(:,3);
    end
end

%% Assemble coefficients Matrix 
K = W-V*tan(phi);
%% Solve the Vortex Strengths
gamma = K\(-Uinf*alpha*ones(Num_panels,1));

% Order Vortex Strengths in a matrix (for the Startboard wing)
vortex=0;
for i=1:npx
    for j=1:npy  
       vortex=vortex+1;
       V_strength(i,j) = gamma(vortex);
           
    end
         
end

%% Lift Coefficient

CL  = (4*s*sum(V_strength(:)))....
    /(S*Uinf*npy);

%% Spanwise Lift Distribution
% ----------------------------------------------------------------------- %
syms y
% Linear chord symbolic function (m)


for j =1:npy % Solve only the Starboard wing due to symmetry  

   y(j)= Panels(1,j).ym; % Local y coordinate in adimensional form    
   c(j)=0.5* (dx_A(j)+dx_B(j))*npx;  % Local chord at the coordinate y (linear chord function)
   cl_nth(j)= 2*sum((V_strength(:,j))/(Uinf*c(j))); % Columwise sumatory of CL_nth  
end
%plotting
figure(2)
plot(y,cl_nth,'o-b')
xlabel(' Local y coordinate')
ylabel('Lift Coefficient')
hold on
plot (-y,cl_nth,'o-b')

count=1;
for al=1:10
alphaa=al;
Gamma= K\(-Uinf*alphaa*ones(Num_panels,1));
vortex=0;
 for i=1:npx
    for j=1:npy  
       vortex=vortex+1;
       Vstrength(i,j) = Gamma(vortex);
           
    end
         
 end

CL_al=(4*s*sum(Vstrength(:)))/(S*Uinf*npy);
cl_plot(count)=CL_al;
count=count+1;
end 
%plotting cl vs alpha
figure(3)
alpha_plot=[1:1:10];
plot(alpha_plot,cl_plot,'-o')
grid on
title('Wing lift coefficient (CL) vs Angle of Attack (alpha)')
ylabel('Wing Lift Coefficient, CL')
xlabel('Angle of Attack, alpha [deg]')

 

