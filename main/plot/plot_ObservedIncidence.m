function [] = plot_ObservedIncidence()
%Summary of this function goes here

  I = getI();
  days = 120+(1:7:43*7+1);
  x = days;
  y = I;
  %figure;
  %plot(x,y/(3.6*100000));
  %plot(x,y/(3.6*10));
  %xlim=120+[0:306];
  %hold on;
  %B.FaceColorData = uint8(255*[1;0;0;0.3]); 
  %alpha(get(B,'children'),0.2);
  %set(ch,'facea',.1)
[AX,H1,H2] = plotyy(x,-1+zeros(1,length(x)),x,y', 'plot');
set(H2,'linewidth',2);
set(AX(1),'xlim',[120,120+30.5*10-2]);
set(AX(2),'xlim',[120,120+30.5*10-2]);
set(AX(1),'ylim',[0,120]);
set(AX(1),'YTick',[0:30:120]);
set(AX(1),'ycolor','b'); 
%set(get(AX(1),'Ylabel'),'String','Susceptible, cumulative incidence (%)', 'FontSize', 12)  
set(get(AX(2),'Ylabel'),'String','Labotory confirmed cases (weekly)', 'FontSize', 12)  
function i = getI()
i = [
2
0
1
9
47
23
72
275
384
239
343
558
930
1261
1344
1082
1630
1760
1659
2400
2914
3326
1660
1013
469
257
234
123
115
142
296
319
260
249
237
258
210
214
202
150
113
98
81
64
];
end
end
