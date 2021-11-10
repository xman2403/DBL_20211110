% use linear regression to determine the slope of linear region above SWI 
% for CCR/LH % 12/03/2021
% read CCR_arranfed_2801_20210311.m
% read Bulk and penetration 
% try to determine SWI by curve fitting
% the number of profile at that day
 row = 1;  
 DBL = zeros (1,74);
 y_DBL_effective = zeros (1,74);
 sed = zeros (50,74);
 C_swi = zeros (1,74); % cswi % 16/02/2021
 water = zeros (50,74);
 Untitled = DO_arranged';
 height_arranged = height_arranged';
 height_arranged(:,115:130)=NaN;
 height_arranged (:,:) = 99 - height_arranged (:,:);
 
 for row = 22:46
   
 % Bulk is the x axis where bulk water ends
   Bulk = bulk(1,row);
   penetration = Penetration (1,row);
  
  % below is dx/dDO;
    [curve,s] = polyfit(Untitled(row,(Bulk):(penetration)),height_arranged(row,(Bulk):(penetration)),6);
  yy2 = polyval(curve,Untitled(row,(Bulk):(penetration)));
 % yy2_b_p = yy2 ((Bulk):(penetration));
  
  y = height_arranged(row,(Bulk):(penetration));
  
      yresid = y - yy2;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1)*var(y);
    R_square = 1 - SSresid/SStotal;
  
   % below R2 only for fitting bulk to penetration
  % R_square = 1 - (s.normr/norm(y - mean(y)))^2;
  
    xx = Untitled(row,(Bulk):(penetration));  
  % d1
  %  diff1 = 4*curve(1).*xx.^3 + 3*curve(2).*xx.^2 + 2*curve(3).*xx + curve(4);
  % diff1 = 5*curve(1).*xx.^4 + 4*curve(2).*xx.^3 + 3*curve(3).*xx.^2 + 2*curve(4).*xx+curve(5);
    diff1 = 6*curve(1).*xx.^5 + 5*curve(2).*xx.^4 + 4*curve(3).*xx.^3 + 3*curve(4).*xx.^2 + 2*curve(5).*xx+curve(6);
  % diff1 = 8*curve(1).*xx.^7 + 7*curve(2).*xx.^6 + 6*curve(3).*xx.^5 + 5*curve(4).*xx.^4 + 4*curve(5).*xx.^3 + 3*curve(6).*xx.^2 + 2*curve(7).*xx+curve(8);
  % d2 
  % diff2 = 7*8*curve(1).*xx.^6 + 6*7*curve(2).*xx.^5 + 5*6*curve(3).*xx.^4 + 4*5*curve(4).*xx.^3 + 3*4*curve(5).*xx.^2 + 2*3*curve(6).*xx.^1 + 2*curve(7);
% 8-order poly
 %  diff2 = 3*4*curve(1).*xx.^2 + 2*3*curve(2).*xx + 2*curve(3);
 %  diff2 = 4*5*curve(1).*xx.^3 + 3*4*curve(2).*xx.^2 + 2*3*curve(3).*xx + 2*curve(4);
    diff2 = 6*5*curve(1).*xx.^4 + 4*5*curve(2).*xx.^3 + 3*4*curve(3).*xx.^2 + 2*3*curve(4).*xx + 2*curve(5);
 
 %% determine SWI % 11/12/2020
   j = 1;
   col = 1;
   for j = (Bulk):(penetration-1)
      if    (diff2(j-Bulk+1)*diff2(j-Bulk+2)<0)
         DBL(col,row) = height_arranged(row,j);
         col = col+1;
      end
   end
   
 DBL(DBL==0)=NaN;
   
 [row3,SWI] = find(height_arranged(row,:) ==max(DBL(:,row)));
 [row3,DBL_upper] = find(height_arranged(row,:) ==min(DBL(:,row)));
 
  sed (row,1:30) = Untitled (row, (SWI+1):(SWI+30));
  C_swi(1,row) = Untitled (row,SWI-1);
  water (row,1:SWI) = Untitled (row,1:SWI);
 
  %  figure; plot (height_arranged (row,(1):(penetration)),Untitled(row,(1):(penetration)),'r--o'); hold on; 
 % plot(yy2,Untitled (row,(Bulk):(penetration)),'b'); hold on;
 % plot([max(DBL(:,row)),max(DBL(:,row))],[Untitled(row,Bulk),Untitled(row,penetration)]); hold on;
 % plot([min(DBL(:,row)),min(DBL(:,row))],[Untitled(row,Bulk),Untitled(row,penetration)]); xlim([height_arranged(row,Bulk), height_arranged(row,penetration)]);
 
 [p,q] = polyfit (height_arranged(row,DBL_upper:SWI),Untitled(row,DBL_upper:(SWI)),1);
    yy = polyval(p,height_arranged(row,DBL_upper:SWI));
  %  line([height_arranged(row,1),height_arranged(row,SWI)],[yy(1,(SWI-DBL_upper+1))-p(1,1)*(height_arranged (row,SWI)-height_arranged (row,1)),yy(1,SWI-DBL_upper+1)],'color','k');  hold on;
     
     % read Cbulk from file % 01/02/2021
     C_bulk = mean(Untitled (row,1:Bulk));
     y_DBL_effective(1,row) = (yy(1,(SWI-DBL_upper+1)) -C_bulk)/p(1,1);
     
     % draw Cbulk and DBL upper boundary (effective) % 02/02/2021
  %   line([height_arranged(row,2),height_arranged(row,SWI)],[C_bulk,C_bulk],'color','c'); xlim([86,105]); hold on;
  %   line([(max(DBL(:,row))-y_DBL_effective(1,row)),(max(DBL(:,row))-y_DBL_effective(1,row))],[0,10]); ylim([0,12])
  %   view(90,90);  % convert x and y axis % 08/03/2021
     
  row2 = row+200; % LH has more than 100 profiles for 1624 % 23/03/2021
%  fig = gcf;
%  fig.PaperUnits = 'inches';
%  fig.PaperPosition = [0 0 3 3];
%  text (97,4,num2str(R_square));
%  print(num2str(row2),'-dpng','-r0');
%  saveas(gcf,[num2str(row2),'.fig']);
%  close(figure(gcf));
   
   
%  figure; 
%  plot (height_arranged (row,(Bulk):(penetration)), diff1,'m--o');hold on;
%  plot (height_arranged (row,(Bulk):(penetration)),diff2,'k-.^'); xlim([86,105]);
%  fig = gcf;
%  fig.PaperUnits = 'inches';
%  fig.PaperPosition = [0 0 3 3];
%  print(num2str(row),'-dpng','-r0');
%  saveas(gcf,[num2str(row),'.fig']);
%  close(figure(gcf));
 end % row cycle