%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA107
% Project Title: Implementation of Differential Evolution (DE) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function z=ESinv(x)
global data PL Cost_Op 
  
%% Example 3A Allen J. Wood
%   Fp=   [data(1,1).*x(1,1).^2+data(1,2)*x(1,1)+data(1,3);                        
%          data(2,1).*x(1,2).^2+data(2,2)*x(1,2)+data(2,3);
%          data(3,1).*x(1,3).^2+data(3,2)*x(1,3)+data(3,3);];
%   z=sum(Fp);
%   lam = abs(sum(x)-PL);
%   z=sum(Fp)+10.*lam; %Pgen fitness function

%% Operation Cost - Yearly
No_BattCost = evalin('base','Cost_Op');  
Konv = evalin('base','konvergen');  
if Konv == 1 
No_BattCostYear = 365*No_BattCost*Konv;  %yearly
else
"Not convergent"
No_BattCostYear = 9.99*10^10;
end
%% Li-Ion Equivalent Uniform Annual Cost
BattInv = sum(0.209*72*x(1,1:4)*0.1295);

z = BattInv + No_BattCostYear;
    
end

