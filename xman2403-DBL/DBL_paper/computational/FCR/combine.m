% use linear regression to determine the slope of linear region above SWI 
% 09/11/2020
% read water_201005_201027.csv;
% read 201111.csv % read the whole profile for the epilimnion core
% try to determine SWI by curve fitting
% the number of profile at that day
 row = 2;  
 DBL = zeros (1,50);
 y_DBL_effective = zeros (1,50);
 sed = zeros (50,50);
 C_swi = zeros (1,50); % cswi % 16/02/2021
 water = zeros (50,50);
 Untitled = epi201206;
 
 for row = 2:38
   
 % Bulk is the x axis where bulk water ends
   Bulk = bulk_20210125 (1,row);
  
  penetration = Penetration (1,row);
  
  % below is dx/dDO;
    [curve,s] = polyfit(Untitled(row,(Bulk):(penetration)),Untitled(1,(Bulk):(penetration)),6);
  yy2 = polyval(curve,Untitled(row,(Bulk):(penetration)));
 % yy2_b_p = yy2 ((Bulk):(penetration));
  
  y = Untitled(1,(Bulk):(penetration));
  
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
         DBL(col,row) = Untitled(1,j);
         col = col+1;
      end
   end
   
 DBL(DBL==0)=NaN;
   
 [row3,SWI] = find(Untitled ==max(DBL(:,row)));
 [row3,DBL_upper] = find(Untitled ==min(DBL(:,row)));
 
  sed (row,1:50) = Untitled (row, (SWI+1):(SWI+50));
  C_swi(1,row) = Untitled (row,SWI);
  water (row,1:SWI) = Untitled (row,1:SWI);
 
%    figure; plot (Untitled (1,(1):(penetration)),Untitled(row,(1):(penetration)),'r--o'); hold on; 
%  plot(yy2,Untitled (row,(Bulk):(penetration)),'b'); hold on;
%  plot([max(DBL(:,row)),max(DBL(:,row))],[Untitled(row,Bulk),Untitled(row,penetration)]); hold on;
%  plot([min(DBL(:,row)),min(DBL(:,row))],[Untitled(row,Bulk),Untitled(row,penetration)]); xlim([Untitled(1,Bulk), Untitled(1,penetration)]);
 
 [p,q] = polyfit (Untitled(1,DBL_upper:SWI),Untitled(row,DBL_upper:(SWI)),1);
    yy = polyval(p,Untitled(1,DBL_upper:SWI));
 %   line([Untitled(1,1),Untitled(1,SWI)],[yy(1,(SWI-DBL_upper+1))-p(1,1)*(Untitled (1,SWI)-Untitled (1,1)),yy(1,SWI-DBL_upper+1)],'color','k');  hold on;
     
     % read Cbulk from file % 01/02/2021
     
     y_DBL_effective(1,row) = (yy(1,(SWI-DBL_upper+1)) -C_bulk(1,row-1))/p(1,1);
     
     % draw Cbulk and DBL upper boundary (effective) % 02/02/2021
  %   line([Untitled(1,2),Untitled(1,SWI)],[C_bulk(1,row-1),C_bulk(1,row-1)],'color','c'); xlim([94,100]); hold on;
  %   line([(max(DBL(:,row))-y_DBL_effective(1,row)),(max(DBL(:,row))-y_DBL_effective(1,row))],[0,10]); ylim([0,12])
   %  view(90,90);  % convert x and y axis % 08/03/2021
     
  row2 = row+100;
%  fig = gcf;
%  fig.PaperUnits = 'inches';
%  fig.PaperPosition = [0 0 3 3];
%  text (97,4,num2str(R_square));
%  print(num2str(row2),'-dpng','-r0');
%  saveas(gcf,[num2str(row2),'.fig']);
%  close(figure(gcf));
   
   
%  figure; 
%  plot (Untitled (1,(Bulk):(penetration)), diff1,'m--o');hold on;
%  plot (Untitled (1,(Bulk):(penetration)),diff2,'k-.^'); xlim([94,100]);
%  fig = gcf;
%  fig.PaperUnits = 'inches';
%  fig.PaperPosition = [0 0 3 3];
%  print(num2str(row),'-dpng','-r0');
 % saveas(gcf,[num2str(row),'.fig']);
 % close(figure(gcf));
 end % row cycle