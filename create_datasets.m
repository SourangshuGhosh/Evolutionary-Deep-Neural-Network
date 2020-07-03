function create_datasets(problem,num_points)
% problem = 'DTLZ4';
% num_points = 100;
switch problem
	case 'ZDT1'
		x = lhsdesign(num_points,30);
		f1 = x(:,1);
		g2_1 = x(:,2:30);
		g2_2 = sum(g2_1,2);
		%g2 = 1+((9*g2_2)/29);
        g2=1;
		f2 = g2 .* (1- sqrt((f1./g2)));
		xlswrite([problem '_' num2str(num_points)],[x,f1,f2],1,'A1');
	case 'ZDT2'
		x = lhsdesign(num_points,30);
		f1 = x(:,1);
		g2_1 = x(:,2:30);
		g2_2 = sum(g2_1,2);
		%g2 = 1+((9*g2_2)/29);
        g2=1;
		f2 = g2 .* (1- power((f1./g2),2));
		xlswrite([problem '_' num2str(num_points)],[x,f1,f2],1,'A1');
	case 'ZDT3'
		%x = lhsdesign(num_points,30);
        x = rand(num_points,30);
        f1 = x(:,1);
		g2_1 = x(:,2:30);
		g2_2 = sum(g2_1,2);
		%g2 = 1+((9*g2_2)/29);
        g2=1;
		f2 = g2 .* (1 - (sqrt((f1 ./ g2)) - ((f1 ./ g2) .* (sin(10*pi*f1)))));
		xlswrite([problem '_' num2str(num_points)],[x,f1,f2],1,'A1');
	case 'ZDT4'
		x = lhsdesign(num_points,10);
		f1 = x(:,1);
		g2_1 = x(:,2:10);
		g2_1 = power(g2_1,2) - 10*cos(4*pi*g2_1);
		g2_2 = sum(g2_1,2);
		%g2 = 1+ 90 + g2_2;
        g2=1;
		f2 = g2 .* (1-(sqrt(f1./g2)));
		xlswrite([problem '_' num2str(num_points)],[x,f1,f2],1,'A1');
	case 'ZDT6'
		x = lhsdesign(num_points,10);
		g1 = x(:,1);
		f1 = 1 - exp(-4*g1) .* power(sin(6*pi*g1),6);
		g2_1 = x(:,2:10);
		g2_2 = sum(g2_1,2);
		%g2 = 1 + 9 * power((g2_2./9),0.25);
        g2=1;
		f2 = g2 .* (1 - power((f1./g2),2));
		xlswrite([problem '_' num2str(num_points)],[x,f1,f2],1,'A1');
	case 'DTLZ1'
		k = 5;
		M = 3;
		n = (M-1) + k; %this is the default
		x = lhsdesign(num_points,n);
		xm = x(:,n-k+1:end); %xm contains the last k variables
		g = 100*(k + sum((xm - 0.5).^2 - cos(20*pi*(xm - 0.5)),2));
		f(:,1) = 1/2*prod(x(:,1:M-1),2).*(1 + g);
		for ii = 2:M-1
			f(:,ii) = 1/2*prod(x(:,1:M-ii),2).*(1 - x(:,M-ii+1)).*(1 + g);
		end
		f(:,M) = 1/2*(1 - x(:,M)).*(1 + g);
		xlswrite([problem '_' num2str(num_points)],[x,f],1,'A1');
	case 'DTLZ2'
		k = 10;
		M = 3;
		n = (M-1) + k; %this is the default
		x = lhsdesign(num_points,n);
		xm = x(:,n-k+1:end); %xm contains the last k variables
		g = sum((xm - 0.5).^2, 2);
		f(:,1) = (1 + g).*prod(cos(pi/2*x(:,1:M-1)),2);
		for ii = 2:M-1
			f(:,ii) = (1 + g) .* prod(cos(pi/2*x(:,1:M-ii)),2) .* ...
				sin(pi/2*x(:,M-ii+1));
		end
		f(:,M) = (1 + g).*sin(pi/2*x(:,1));
		xlswrite([problem '_' num2str(num_points)],[x,f],1,'A1');
	case 'DTLZ3'
		k = 10;
		M = 3;
		n = (M-1) + k; %this is the default
		x = lhsdesign(num_points,n);
		xm = x(:,n-k+1:end); %xm contains the last k variables
		g = 100*(k + sum((xm - 0.5).^2 - cos(20*pi*(xm - 0.5)),2));
		f(:,1) = (1 + g).*prod(cos(pi/2*x(:,1:M-1)),2);
		for ii = 2:M-1
			f(:,ii) = (1 + g) .* prod(cos(pi/2*x(:,1:M-ii)),2) .* ...
				sin(pi/2*x(:,M-ii+1));
		end
		f(:,M) = (1 + g).*sin(pi/2*x(:,1));
		xlswrite([problem '_' num2str(num_points)],[x,f],1,'A1');
	case 'DTLZ4'
		k = 10;
		M = 3;
		n = (M-1) + k; %this is the default
		alpha = 100; %suggested value is 100
		x = lhsdesign(num_points,n);
		xm = x(:,n-k+1:end); %xm contains the last k variables
		g = sum((xm - 0.5).^2, 2);
		f(:,1) = (1 + g).*prod(cos(pi/2*x(:,1:M-1).^alpha),2);
		for ii = 2:M-1
			f(:,ii) = (1 + g) .* prod(cos(pi/2*x(:,1:M-ii).^alpha),2) .* ...
				sin(pi/2*x(:,M-ii+1).^alpha);
		end
		f(:,M) = (1 + g).*sin(pi/2*x(:,1).^alpha);
		xlswrite([problem '_' num2str(num_points)],[x,f],1,'A1');
	case 'DTLZ5'
		k = 10;
		M = 3;
		n = (M-1) + k; %this is the default
		x = lhsdesign(num_points,n);
		xm = x(:,n-k+1:end); %xm contains the last k variables
		g = sum((xm - 0.5).^2, 2);	
		theta(1,:) = pi/2*x(:,1);
		gr = g(:,ones(M-2,1)); %replicates gr for the multiplication below
		theta(:,2:M-1) = pi./(4*(1+gr)) .* (1 + 2*gr.*x(:,2:M-1));
		f(:,1) =(1 + g).*prod(cos(theta(:,1:M-1)),2);
		for ii = 2:M-1
			f(:,ii) = (1 + g) .* prod(cos(theta(:,1:M-ii)),2) .* ...
				sin(theta(:,M-ii+1));
		end
		f(:,M) = (1 + g).*sin(theta(:,1));
		xlswrite([problem '_' num2str(num_points)],[x,f],1,'A1');
	case 'DTLZ6'
		k = 10;
		M = 3;
		n = (M-1) + k; %this is the default
		x = lhsdesign(num_points,n);
		xm = x(:,n-k+1:end); %xm contains the last k variables
		g = sum(xm.^0.1, 2);	
		theta(1,:) = pi/2*x(:,1);
		gr = g(:,ones(M-2,1)); %replicates gr for the multiplication below
		theta(:,2:M-1) = pi./(4*(1+gr)) .* (1 + 2*gr.*x(:,2:M-1));
		f(:,1) =(1 + g).*prod(cos(theta(:,1:M-1)),2);
		for ii = 2:M-1
			f(:,ii) = (1 + g) .* prod(cos(theta(:,1:M-ii)),2) .* ...
				sin(theta(:,M-ii+1));
		end
		f(:,M) = (1 + g).*sin(theta(:,1));
		xlswrite([problem '_' num2str(num_points)],[x,f],1,'A1');
	case 'DTLZ7'
		k = 20;
		M = 3;
		n = (M-1) + k; %this is the default
		x = lhsdesign(num_points,n);
		xm = x(:,n-k+1:end); %xm contains the last k variables
		g = 1 + 9/k*sum(xm,2);
		f(:,1:M-1) = x(:,1:M-1);
		% The last function requires another auxiliar variable
		gaux = g(:,ones(M-1,1)); %replicates the g function
		h = M - sum(f./(1+gaux).*(1 + sin(3*pi*f)),2);
		f(:,M) = (1 + g).*h;
		xlswrite([problem '_' num2str(num_points)],[x,f],1,'A1');

end