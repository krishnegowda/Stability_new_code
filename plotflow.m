function plotflow(y, flow)
  yy = y(2:end-1);  % the flow doesn't contain the boundaries
  n = length(yy);
  v = flow(1:n);
  vort = flow(n+1:end);
  plot(yy, real(v), 'k');
  hold on
  plot(yy, imag(v), 'k--');
  plot(yy, real(vort), 'b');
  plot(yy, imag(vort), 'b--');
  hold off