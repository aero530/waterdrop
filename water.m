%Implementation of the 2D wave equation
%                     ©Julius Piso 2014
% http://www.piso.at/julius/index.php/projects/programming/13-2d-wave-equation-in-octave

s = 128; %array size (spacial resolution of the simulation)

%create arrays
p = zeros(s,s);  %past
c = zeros(s,s);  %current
f = zeros(s,s);  %future

speed = 10;  			%propagation speed
dt = 0.0001; 			%timestep (small!!) (time resolution of the simulation)
dx = 0.01; 				%distance between elements
r = speed * dt / dx; 
stopTime = 10;        	%stop time

%for dynamic boundary
n = 300;

%initial loop
for i = s/2 - s/16 : s/2 + s/16
	%initial conditions
	%c(i,1:s) = 0.5 * sin(2*pi/s*[1:s]) * sin(i *2*pi/s);
	c(i, s/2 - s/16 : s/2 + s/16) = - 2 * cos(0.5 * 2 * pi / (s/8) * [s/2 - s/16 : s/2 + s/16]) * cos(0.5 * 2 * pi / (s/8) * i);
	p(i,1:s) = c(i,1:s);
end

%main loop
for t = 0 : dt : stopTime
	%wave equation
	f(2 : s-1,2:s-1) = 2 * c(2 : s-1,2:s-1) - p(2 : s-1,2:s-1) + r ^ 2 * ( c(1 : s-2,2:s-1) + c(3 : s,2:s-1) + c(2 : s-1,1:s-2) + c(2 : s-1,3:s) - 4 * c(2 : s-1,2:s-1) );
	
	%body force (dynamic boundary)
	%f(s/2+s/4-1 : s/2+s/4+1 , 1:2) = 1.5 * sin(t*n);
	%f(s/2-s/4-1 : s/2-s/4+1 , 1:2) = 1.5 * sin(t*n);
	%f(2 : s-1,1:2) = 1;
	
	%transparent boundaries (no reflection)
	f(2:s-1,1) = (2 * c(2:s-1,1) + (r-1) * p(2:s-1,1) + 2*r*r*(c(2:s-1,2) + c(3:s,1) + c(1:s-2,1) - 3 * c(2:s-1,1)))/(1+r);   % Y:1
	f(2:s-1,s) = (2 * c(2:s-1,s) + (r-1) * p(2:s-1,s) + 2*r*r*(c(2:s-1,s-1) + c(3:s,s) + c(1:s-2,s) - 3 * c(2:s-1,s)))/(1+r); % Y:s
	f(1,2:s-1) = (2 * c(1,2:s-1) + (r-1) * p(1,2:s-1) + 2*r*r*(c(2,2:s-1) + c(1,3:s) + c(1,1:s-2) - 3 * c(1,2:s-1)))/(1+r);   % X:1
	f(s,2:s-1) = (2 * c(s,2:s-1) + (r-1) * p(s,2:s-1) + 2*r*r*(c(s-1,2:s-1) + c(s,3:s) + c(s,1:s-2) - 3 * c(s,2:s-1)))/(1+r); % Y:s
	
	
	%special corners (transparent boundary)
	f(1,1) = ( 2 * c(1,1) + ( r-1 ) * p(1,1) + 2*r*r* ( c(2,1) + c(1,2) - 2*c(1,1) )) / (1+r);  % X:1 ; Y:1
    f(s,s) = ( 2 * c(s,s) + (r-1) * p(s,s) + 2*r*r*(c(s-1,s) + c(s,s-1) - 2*c(s,s))) / (1+r);   % X:s ; Y:s
    f(1,s) = ( 2 * c(1,s) + (r-1) * p(1,s) + 2*r*r* ( c(2,s) + c(1,s-1) - 2*c(1,s))) / (1+r);   % X:1 ; Y:s
    f(s,1) = ( 2 * c(s,1) + (r-1) * p(s,1) + 2*r*r* ( c(s-1,1) + c(s,2) - 2*c(s,1))) / (1+r);   % X:s ; Y:1
  
	%update values
	p(1 : s,1:s) = c(1 : s,1:s);
	c(1 : s,1:s) = f(1 : s,1:s);
	
	%plot
	if mod( t / dt , 10 ) == 0
		b = surf(c); 				%surface plot
		%set(b, 'edgecolor','none') %turn mesh edges off
		axis([0 s 0 s -2 2])  		%axis scale
		caxis([-1 1])         		%color axis scale
		pause ( .0001 )
	end
	
end
