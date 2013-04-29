rmin = 1;
rmax = 4;
rstep = 0.01;

figure
hold on
for r = rmin:rstep:rmax
   
   xinit = 0.01;
   
   for j = 1:50
      xnext = r * xinit * (1 - xinit);
      xinit = xnext;
   end
   
   for j = 1:16
      xnext = r * xinit * (1 - xinit);
      xinit = xnext;
      scatter(r, xnext, '.b', 'SizeData', 2);
   end

end   
hold off