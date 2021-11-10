% calculate DBL thickness using k1 and Cswi % 03/11/2020
% load 'water_201103.mat';
  Cb = C_bulk(1:41);  % (mg/L)
  Cswi = C_swi(1:41);
  k = 223/86400;   % m1 for FCR epi core (divide 86400 to get s-1)
 % k = k_m4_201103./86400;
  Dw = 1.34e-9;    % m2/s
  Ds = 1.27e-9;
% DBL_sim = 1000*(Cswi./Cb-1).*Dw./((k.*Ds).^0.5);   % mm
  DBL_sim = 1000*(Cb./Cswi-1).*Dw./((k*Ds)^0.5);   % mm
  
% calculate JO2 at water side % 04/11/2020
JO2_sim = (Dw./DBL_sim).*(Cb-Cswi)*86400*1000000/32;  % mmolm-2d-1