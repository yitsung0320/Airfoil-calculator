% Airfoil Data Calculator
% Develped by Yi Tsung Lee
% First edit 09062019

clc
clear
% ====== File Reading ============== 

%{
Data = importdata('file1.txt',' ',1) %Creating Data "object"
%Data.data is the valuew inside
%Data.textdata is the header now


[x_index,y_index] = size(Data.data);
data = Data.data ;
point_index = x_index ;
%}


%
% Open the file and store the pointer to the file in fid
% 
fid = fopen('file2.txt','rt');
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

[column_index, row_index] = size(data);
point_index = column_index ;
for i = 1:column_index
    xp(1,i)     = data(i,2); % xpoint coordinate
    yp(1,i)     = data(i,3); % ypoint coordinate
    v_vinf(1,i) = data(i,4); % v/vinf data
    
    cf(1,i)     = data(i,7); % local skin friction coefficient
end

% ====== Plot of Airfoil Shape ===============
subplot (2,1,1)
plot(xp(1,:),yp(1,:),'k-')
title('Airfoil SHape ')
xlabel('x/c')
ylabel('y')

% ====== cp calculation ============
% cp = 1-(v/vinf)^2

for i = 1:point_index

    cp(1,i) = 1-v_vinf(1,i).^2;

end


% ====== Find Stagnation Point =======
% at stagnation pt : cp = 1; using iteration to find the closet cp point
% the Cp examination is not in use thiis time.
%{
cp0 = 1-cp; % examine the 1-cp to see if it is close to 0
stag_pt_index =find(cp0 == min(cp0));
%}

% another way to preciselyfind the stagnation pt is to find where V is change the
% sign
for i = 1: point_index-1
    
  if v_vinf(1,i)*v_vinf(1,i+1) <= 0  
    stag_pt = i+(0-v_vinf(1,i)/(v_vinf(1,i+1)-v_vinf(1,i)));  % this is the precise stag pt : a decimal
    stag_pt_index2 = i;  % the index of pt right before before stag
    stag_ratio = stag_pt-stag_pt_index2;  % the position in % of the stag in the panel
    break %break of for
  end
  
end

% ======= interpolating data froom stagnation pt ============

xp_stag = xp(1,stag_pt_index2) + stag_ratio*(xp(1,stag_pt_index2+1)-xp(1,stag_pt_index2));
yp_stag = yp(1,stag_pt_index2) + stag_ratio*(yp(1,stag_pt_index2+1)-yp(1,stag_pt_index2));
cp_stag = cp(1,stag_pt_index2) + stag_ratio*(cp(1,stag_pt_index2+1)-cp(1,stag_pt_index2));
cf_stag = cf(1,stag_pt_index2) + stag_ratio*(cf(1,stag_pt_index2+1)-cf(1,stag_pt_index2));

% ====== Force from pressure ==================

% == Fundamental Assumption == : 
v_inf = 100;    % Noted the ration of v_inf and p_inf has to be 
                % large enough so that the decimal error frpm cp will be
                % calculated
p_inf = 1;
pho = 1;

for i =1:point_index-1
    
    cp_ave(1,i) = (cp(1,i+1) + cp(1,i))/2 ;
    px(1,i) = -1*(p_inf +1/2*pho*v_inf^2*cp_ave(1,i))*(yp(1,i+1)-yp(1,i));
    py(1,i) =    (p_inf +1/2*pho*v_inf^2*cp_ave(1,i))*(xp(1,i+1)-xp(1,i));

end
X_p = sum(px); % X force from pressure
Y_p = sum(py); % X force from Pressure

Cx_p = X_p/(1/2*pho*v_inf^2)
Cy_p = Y_p/(1/2*pho*v_inf^2)
% ===== Force from Friction ===================

for i =1:point_index-1
   
    cf_ave(1,i) = (cf(1,i+1) + cf(1,i))/2 ;
    if i < stag_pt_index2
      fx(1,i) = -1/2*pho*v_inf^2*cf_ave(1,i)*(xp(1,i+1)-xp(1,i));
      fy(1,i) = -1/2*pho*v_inf^2*cf_ave(1,i)*(yp(1,i+1)-yp(1,i));
    else
      fx(1,i) =  1/2*pho*v_inf^2*cf_ave(1,i)*(xp(1,i+1)-xp(1,i));
      fy(1,i) =  1/2*pho*v_inf^2*cf_ave(1,i)*(yp(1,i+1)-yp(1,i));
    end  %end of if
    
end
 
X_f = sum(fx)...
      -1/2*pho*v_inf^2*(cf(1,stag_pt_index2)+cf_stag)/2*(xp_stag-xp(1,stag_pt_index2))...
      +1/2*pho*v_inf^2*(cf(1,stag_pt_index2+1)+cf_stag)/2*(xp(1,stag_pt_index2+1)-xp_stag); % X force from friction
Y_f = sum(fy)...
      -1/2*pho*v_inf^2*(cf(1,stag_pt_index2)+cf_stag)/2*(yp_stag-yp(1,stag_pt_index2)) 
      +1/2*pho*v_inf^2*(cf(1,stag_pt_index2+1)+cf_stag)/2*(yp(1,stag_pt_index2+1)-yp_stag); % Y force from friction
Cx_f = X_f/(1/2*pho*v_inf^2)
Cy_f = Y_f/(1/2*pho*v_inf^2)

% ===== SUmmary of Force =====================
Cx = (X_p + X_f)/(1/2*pho*v_inf^2*1)  % Force Coeff in x
Cy = (Y_p + Y_f)/(1/2*pho*v_inf^2*1)  % Force Coeff in y

%{
prompt = 'what Angle of Attack is the airfoil in (degree)?';
alpha = input(prompt)/180*pi ;
%}

alpha = 8/180*pi ;

Cl = Cy*cos(alpha) - Cx*sin(alpha)  % Airfoil Lift coeff
Cd = Cy*sin(alpha) + Cx*cos(alpha)  % Airfoil Drag coeff
 
Cl_p = Cy_p*cos(alpha) - Cx_p*sin(alpha) % pressure contribution of Cl
Cd_p = Cy_p*sin(alpha) + Cx_p*cos(alpha) % pressure contribution of cd

Cl_f = Cy_f*cos(alpha) - Cx_f*sin(alpha) % friction contribution of Cl
Cd_f = Cy_f*sin(alpha) + Cx_f*cos(alpha) % friction contribution of Cd

% ======Pitching Moment Coefficient =========
% 1/4 Cchord Cm 

for i = 1: point_index-1
    
     xp_ave(1,i) =(xp(1,i+1)+xp(1,i))/2;
     yp_ave(1,i) =(yp(1,i+1)+yp(1,i))/2; 
     mp(1,i) =  -py(1,i)*(xp_ave(1,i)-1/4) ...
               +px(1,i)*(yp_ave(1,i));
     mf(1,i) = -fy(1,i)*(xp_ave(1,i)-1/4)...
               +fx(1,i)*(yp_ave(1,i));
end

Cm_p = sum(mp)/((1/2*pho*v_inf^2)*1*1)
Cm_f = sum(mf)/((1/2*pho*v_inf^2)*1*1)
Cm = sum(mp+mf)/((1/2*pho*v_inf^2)*1*1)

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