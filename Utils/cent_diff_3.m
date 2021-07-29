function dx = cent_diff_3(x,h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the derivative using the central difference formula
% x is the signal to differentiate
% h is the step size (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(x,1) == 1
  x = x';
end

dx = zeros(size(x));
for k = 2:length(x)-1
    dx(k,:) = (x(k+1,:)-x(k-1,:))./(2*h);
end

dx(1,:) = (x(2,:)-x(1,:))./h;
dx(end,:) = (x(end,:)-x(end-1,:))./h;

