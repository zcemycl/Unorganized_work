function lineObj = animInit(height)
hold on
axis equal
axis([-1.2 0.3 0 0.25])
plot([-1.2 -0.5 -0.5],[height height 0],...
         'k','LineWidth',3);

%link line objects

lineObj.h0T = line(0,0,'color','k','LineWidth',2);
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
lineObj.h5T = line(0,0,'color','k','LineWidth',2);

end 