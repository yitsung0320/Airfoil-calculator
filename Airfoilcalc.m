% MAE 551 HW1 Airfoil Calculator
% Airfoil Data Calculator for hw1 all problem by Yi Tsung Lee
% version eddited on 9/26/2019

clc
clear

% =========== Code stat ==========================
prompt = 'Type the number of the files to read in : ';
file_numbers = input(prompt);

for z = 1 : file_numbers
prompt = 'type the name include the file type to read in :';
filename = input(prompt,'s');  % record as string
if isfile(filename)
else
    filename = input(prompt,'s');
end

% ===== read the file by the filename ============
fid = fopen(filename,'rt');
%
% fgetl(fid) will read the first line of the file
% The information is not used
%
fgetl(fid);
%
% Now read the data into a variable named 'data'
% The following command will read the data into a 7-by-N matrix, 
% where N is the number of rows in the file
% The file will always have 7 columns, but the number of rows
% is unknown
%
data = fscanf(fid,'%f',[7, Inf]);
%
% Transpose the 'data' matrix to make it N-by-7 in size
%
data = data';
%
% Let xte be the x-value of the TE. This is simply the 
% number in the second column, first row
%
xte = data(1,2);
%
% In a data file containing information from a 
%  viscous analysis, the information from the wake 
% (behind the TE) is included, but not useful for our
% program. We need to ignore it.
%
% The following command finds the rows that have 
% information for points in the wake and these rows are
% then removed from 'data'
%
rows_remove = find(data(:,2)>xte);
if (~isempty(rows_remove)),
  data = data(1:rows_remove(1)-1,:);
end
%
% Now this data is ready to use.
%
 
% ====== ectract the value from the file data

[column_index, row_index] = size(data);
point_index(z,1) = column_index ; % number of data point
for i = 1:column_index
    xp(z,i)     = data(i,2); % xpoint coordinate
    yp(z,i)     = data(i,3); % ypoint coordinate
    v_vinf(z,i) = data(i,4); % v/vinf data
    
    cf(z,i)     = data(i,7); % local skin friction coefficient
end

% ====== Find Stagnation Point =======

for i = 1: point_index(z,1)-1
    
  if v_vinf(z,i)*v_vinf(z,i+1) <= 0
    
    % exact position of stag pt in decimal  
    stag_pt(z,1) = i + (0 - v_vinf(z,i)/(v_vinf(z,i+1) - v_vinf(z,i)));
    % the index of pt right before before stag
    stag_pt_index(z,1) = i;  
    % the position in % of the stag in the panel
    stag_ratio(z,1) = stag_pt(z,1)-stag_pt_index(z,1);  % the position in % of the stag in the panel
    break %break of for
  end % end of if
  
end 

% ======= interpolating data froom stagnation pt ============

xp_stag(z,1) = xp(z,stag_pt_index(z,1))...
               + stag_ratio(z,1)*(xp(z,stag_pt_index(z,1)+1)-xp(z,stag_pt_index(z,1)));
yp_stag(z,1) = yp(z,stag_pt_index(z,1))...
               + stag_ratio(z,1)*(yp(z,stag_pt_index(z,1)+1)-yp(z,stag_pt_index(z,1)));
v_vinf_stag(z,1) = v_vinf(z,stag_pt_index(z,1))...
               + stag_ratio(z,1)*(v_vinf(z,stag_pt_index(z,1)+1)-v_vinf(z,stag_pt_index(z,1)));
cf_stag(z,1) = cf(z,stag_pt_index(z,1))...
               + stag_ratio(z,1)*(cf(z,stag_pt_index(z,1)+1)-cf(z,stag_pt_index(z,1)));

% from now on the data above will be inserted on the stag pt 
% shifting the data after stag_pt_index+1 for one row
xp(z,stag_pt_index(z,1)+2:point_index(z,1)+1) = xp(z,stag_pt_index(z,1)+1:point_index(z,1));
yp(z,stag_pt_index(z,1)+2:point_index(z,1)+1) = yp(z,stag_pt_index(z,1)+1:point_index(z,1));
v_vinf(z,stag_pt_index(z,1)+2:point_index(z,1)+1) = v_vinf(z,stag_pt_index(z,1)+1:point_index(z,1));
cf(z,stag_pt_index(z,1)+2:point_index(z,1)+1) = cf(z,stag_pt_index(z,1)+1:point_index(z,1));

% insert the interpolating data at stag_pt_index+1
xp(z,stag_pt_index(z,1)+1) = xp_stag(z,1);
yp(z,stag_pt_index(z,1)+1) = yp_stag(z,1);
v_vinf(z,stag_pt_index(z,1)+1) = v_vinf_stag(z,1);
cf(z,stag_pt_index(z,1)+1) = cf(z,1);

point_index(z,1) = point_index(z,1) +1 ; % after inserting, data pt increase
stag_pt_index(z,1) = stag_pt_index(z,1) +1 ; % now its the real index of stag_pt
% Now we have new data set containing the interpolating data in it 

% ====== cp calculation ============
% cp = 1-(v/vinf)^2

for i = 1:point_index(z,1)
    cp(z,i) = 1-v_vinf(z,i).^2;
end

 % ====== Force coefficient from pressure ==================

for i =1:point_index(z,1)-1  % panel number is pointnumber -1 
    
    cp_ave(z,i) = (cp(z,i+1) + cp(z,i))/2 ;
    cxp(z,i) = -1*cp_ave(z,i)*(yp(z,i+1)-yp(z,i));
    cyp(z,i) =    cp_ave(z,i)*(xp(z,i+1)-xp(z,i));

end

% ===== Force Coefficient from Friction ===================

for i =1:point_index(z,1)-1
   
    cf_ave(z,i) = (cf(z,i+1) + cf(z,i))/2 ;
    if i < stag_pt_index(z,1)
      cxf(z,i) = -cf_ave(z,i)*(xp(z,i+1)-xp(z,i));
      cyf(z,i) = -cf_ave(z,i)*(yp(z,i+1)-yp(z,i));
    else
      cxf(z,i) =  cf_ave(z,i)*(xp(z,i+1)-xp(z,i));
      cyf(z,i) =  cf_ave(z,i)*(yp(z,i+1)-yp(z,i));
    end  %end of if
    
end

% ====== Alpha collection ==========
prompt = 'Type in the alpha(in degree) for this file :';
alpha_ra = input(prompt) ;
alpha(z,1) = alpha_ra*pi/180 ;

% ====== Coefficient calculation =====
Cx_p(z,1) = sum(cxp(z,:)) ;  % pressure contribution of Cx
Cy_p(z,1) = sum(cyp(z,:)) ;  % pressure contribution of Cy
 
Cx_f(z,1) = sum(cxf(z,:)) ;  % friction contribution of Cx
Cy_f(z,1) = sum(cyf(z,:)) ; % friction contribution of Cy

Cx(z,1) = Cx_p(z,1) + Cx_f(z,1); 
Cy(z,1) = Cy_p(z,1) + Cy_f(z,1);

Cl(z,1) = Cy(z,1)*cos(alpha(z,1)) - Cx(z,1)*sin(alpha(z,1));  % Airfoil Lift coeff
Cd(z,1) = Cy(z,1)*sin(alpha(z,1)) + Cx(z,1)*cos(alpha(z,1));  % Airfoil Drag coeff
 
Cl_p(z,1) = Cy_p(z,1)*cos(alpha(z,1)) - Cx_p(z,1)*sin(alpha(z,1)); % pressure contribution of Cl
Cd_p(z,1) = Cy_p(z,1)*sin(alpha(z,1)) + Cx_p(z,1)*cos(alpha(z,1)); % pressure contribution of cd

Cl_f(z,1) = Cy_f(z,1)*cos(alpha(z,1)) - Cx_f(z,1)*sin(alpha(z,1)); % friction contribution of Cl
Cd_f(z,1) = Cy_f(z,1)*sin(alpha(z,1)) + Cx_f(z,1)*cos(alpha(z,1)); % friction contribution of Cd

% ======= Pitching Moment calculation
for i = 1: point_index-1
    
     xp_ave(z,i) =(xp(z,i+1)+xp(z,i))/2;
     yp_ave(z,i) =(yp(z,i+1)+yp(z,i))/2; 
     cmp(z,i) = -cyp(z,i)*(xp_ave(z,i)-1/4)+cxp(z,i)*yp_ave(z,i);
     cmf(z,i) = -cyf(z,i)*(xp_ave(z,i)-1/4)+cxf(z,i)*yp_ave(z,i);
end
Cm_p(z,1) = sum(cmp(z,:));
Cm_f(z,1) = sum(cmf(z,:));
Cm(z,1) = Cm_p(z,1) + Cm_f(z,1);


dataset ={"/n",'Total','Presure','Skin-friction';'Cx',Cx(z,1),Cx_p(z,1),Cx_f(z,1);...
          'Cy',Cy(z,1),Cy_p(z,1),Cy_f(z,1);'Cl',Cl(z,1),Cl_p(z,1),Cl_f(z,1);...
          'Cd',Cd(z,1),Cd_p(z,1),Cd_f(z,1);'Cm',Cm(z,1),Cm_p(z,1),Cm_f(z,1)}
end 

prompt = 'Type the probelm # you are solving :';
problem = input(prompt);

% ===== plot and solution for Porblem 1 =========
if problem == 1
% ====== Plot of c_p ===============
subplot (2,1,1)
plot(xp(1,:),cp(1,:),'B-')
title('c_p distribution of the airfoil')
xlabel('x/c')
ylabel('c_p value')
set(gca,'YDir','reverse')

% ====== Plot of Airfoil Shape ===============
subplot (2,1,2)
plot(xp(1,:),yp(1,:),'k-')
title('Airfoil SHape ')
axis equal
xlabel('x/c')
ylabel('y')   
end %end of if

%========== Problem 2 ===============
if problem == 2
    
dataset2 = {' ','0^o','4^o','8^o','12^o','16^o';...
            'Cl',Cl(1,1),Cl(2,1),Cl(3,1),Cl(4,1),Cl(5,1);...
            'Cd',Cd(1,1),Cd(2,1),Cd(3,1),Cd(4,1),Cd(5,1);...
            'Cd_f',Cd_f(1,1),Cd_f(2,1),Cd_f(3,1),Cd_f(4,1),Cd_f(5,1);...
            'Cd_p',Cd_p(1,1),Cd_p(2,1),Cd_p(3,1),Cd_p(4,1),Cd_p(5,1);...
            'Cm',Cm(1,1),Cm(2,1),Cm(3,1),Cm(4,1),Cm(5,1)}   
figure(1)
subplot(3,2,1)
plot(alpha(:,1)/pi*180,Cl(:,1),'B-')
xlabel('alpha(degree)','FontSize',20)
ylabel('Cl','FontSize',20)

subplot(3,2,2)
plot(alpha(:,1)/pi*180,Cd(:,1),'B-')
xlabel('alpha(degree)','FontSize',20)
ylabel('Cd','FontSize',20)

subplot(3,2,3)
plot(alpha(:,1)/pi*180,Cd_f(:,1),'B-')
xlabel('alpha(degree)','FontSize',20)
ylabel('Cd_f','FontSize',20)

subplot(3,2,4)
plot(alpha(:,1)/pi*180,Cd_p(:,1),'B-')
xlabel('alpha(degree)','FontSize',20)
ylabel('Cd_p','FontSize',20)

subplot(3,2,5)
plot(alpha(:,1)/pi*180,Cm(:,1),'B-')
xlabel('alpha','FontSize',20)
ylabel('Cmc/4','FontSize',20)

subplot(3,2,6)
plot(xp(1,:),yp(1,:),'K-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
axis equal
title('Airfoil shape')

figure(2)
subplot(2,2,1)
plot(xp(1,:),abs(v_vinf(1,:)),'B',xp(2,:),abs(v_vinf(2,:)),'R',...
     xp(3,:),abs(v_vinf(3,:)),'G',xp(4,:),abs(v_vinf(4,:)),'Y',xp(5,:),abs(v_vinf(5,:)),'K')
xlabel('x/c','FontSize',20)
ylabel('V/Vinf','FontSize',20)
legend('alpha = 0','alpha = 4','alpha = 8','alpha = 12','alpha = 16')

subplot(2,2,3)
plot(xp(1,:),yp(1,:),'K-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
axis equal
title('Airfoil shape','FontSize',20)

subplot(2,2,2)
plot(xp(1,:),cp(1,:),'B',xp(2,:),cp(2,:),'R',...
     xp(3,:),cp(3,:),'G',xp(4,:),cp(4,:),'Y',xp(5,:),cp(5,:),'K')
xlabel('x/c','FontSize',20)
ylabel('cp','FontSize',20)
axis([0 1 -inf 1])
set(gca,'YDir','reverse')
legend('alpha = 0','alpha = 4','alpha = 8','alpha = 12','alpha = 16')

subplot(2,2,4)
plot(xp(1,:),yp(1,:),'K-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
axis equal
title('Airfoil shape','FontSize',20)

end

% ======= problem 3 ========
if problem  == 3
figure(1)
subplot(2,1,1)
plot(xp(1,:),cp(1,:),'B',xp(2,:),cp(2,:),'R')
xlabel('x/c','FontSize',20)
ylabel('cp','FontSize',20)
set(gca,'YDir','reverse')
legend('b airfoil','c airfoil')

subplot(2,1,2)
plot(xp(1,:),yp(1,:),'B-',xp(2,:),yp(2,:),'R-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
legend('b airfoil','c airfoil')
axis equal
title('Airfoil shape','FontSize',20)

Cm

end

% == problem 4 ========================

if problem == 4
figure(1)
subplot(2,1,1)
plot(xp(1,:),cp(1,:),'B',xp(2,:),cp(2,:),'R')
xlabel('x/c','FontSize',20)
ylabel('c_p','FontSize',20)
set(gca,'YDir','reverse')
legend('D airfoil','E airfoil')

subplot(2,1,2)
plot(xp(1,:),yp(1,:),'B-',xp(2,:),yp(2,:),'R-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
legend('D airfoil','E airfoil')
axis equal
title('Airfoil shape','FontSize',20)


figure(2)
subplot(2,1,1)
plot(xp(1,:),cf(1,:),'B',xp(2,:),cf(2,:),'R')
xlabel('x/c','FontSize',20)
ylabel('c_f','FontSize',20)
legend('D airfoil','E airfoil')

subplot(2,1,2)
plot(xp(1,:),yp(1,:),'B-',xp(2,:),yp(2,:),'R-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
legend('D airfoil','E airfoil')
axis equal
title('Airfoil shape','FontSize',20)    
end

if problem == 5
    figure(1)
subplot(2,1,1)
plot(xp(1,:),cp(1,:),'B',xp(2,:),cp(2,:),'R')
xlabel('x/c','FontSize',20)
ylabel('c_p','FontSize',20)
set(gca,'YDir','reverse')
legend('Re = 3.0*10^5, AOA = 1.092^o','Re = 3.0*10^6, AOA = 0.798^o')

subplot(2,1,2)
plot(xp(1,:),yp(1,:),'B-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
legend('a airfoil')
axis equal
title('Airfoil shape','FontSize',20)


figure(2)
subplot(2,1,1)
plot(xp(1,:),cf(1,:),'B',xp(2,:),cf(2,:),'R')
xlabel('x/c','FontSize',20)
ylabel('c_f','FontSize',20)
legend('Re = 3.0*10^5,AOA = 1.092^o','Re = 3.0*10^6,AOA = 0.798^o')

subplot(2,1,2)
plot(xp(1,:),yp(1,:),'B-')
xlabel('x/c','FontSize',20)
ylabel('y','FontSize',20)
legend('a airfoil')
axis equal
title('Airfoil shape','FontSize',20)
end