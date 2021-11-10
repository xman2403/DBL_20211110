% use linear regression to determine the slope of linear region above SWI 
% the number of profile is 74
% the traditional effective DBL approach 
% read DO profile mat and height mat before running this code
% read Bulk and penetration mat too
 row = 1;  
 DBL = zeros (1,74);
 y_DBL_effective = zeros (1,74);
 sed = zeros (50,74);
 C_swi = zeros (1,74); % cswi % 16/02/2021
 water = zeros (50,74);

% DO arranged is measured DO profiles; height is lander positions
 Untitled = DO_arranged';
 height_arranged = height_arranged';
 height_arranged(:,115:130)= NaN;
 height_arranged (:,:) = 100 - height_arranged (:,:);
 
 for row = 14:15
  
 % Bulk is the x axis where bulk water ends
   Bulk = bulk(1,row);
     
 % determine SWI 
   
 [row3,SWI] = find(height_arranged(row,:) == 100);
 
 % extrapolate the gradient at SWI

    figure; plot (height_arranged (row,(1):(100)),Untitled(row,(1):(100)),'r--o'); hold on; 
  plot([height_arranged(row,SWI),height_arranged(row,SWI)],[Untitled(row,1),Untitled(row,100)]); xlim([86,100]);
 
 [p,q] = polyfit (height_arranged(row,(SWI-2):SWI),Untitled(row,(SWI-2):(SWI)),1);
    yy = polyval(p,height_arranged(row,(SWI-2):SWI));
    line([height_arranged(row,15),height_arranged(row,SWI)],[yy(1,3)-p(1,1)*(height_arranged(row,SWI)-height_arranged(row,15)),yy(1,3)],'color','k');  hold on;
     
     % Calculate concentration is bulk water
     C_bulk = mean(Untitled (row,1:Bulk));

     % DBL thickness
     y_DBL_effective(1,row) = (yy(1,3) -C_bulk)/(p(1,1));
     
     % Plot Cbulk and DBL upper boundary (effective) 
     line([height_arranged(row,2),height_arranged(row,SWI)],[C_bulk,C_bulk],'color','c'); xlim([86,105]); hold on;
     line([(height_arranged(row,SWI)-y_DBL_effective(1,row)),((height_arranged(row,SWI))-y_DBL_effective(1,row))],[0,10]); ylim([0,12])
     view(90,90);  % convert x and y axis %
     
  % Automatically plot d1d2 and DO profiles, DBL bounds

  row2 = row+100;
  fig = gcf;
  fig.PaperUnits = 'inches';
  fig.PaperPosition = [0 0 3 3];
  print(num2str(row2),'-dpng','-r0');
  saveas(gcf,[num2str(row2),'.fig']);
  close(figure(gcf));
 
 end % row cycle