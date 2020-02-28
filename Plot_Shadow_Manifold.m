% Plot_Shadow_Manifold, Thomas Merkh, tmerkh@ucla.edu, July 25th 2017

function [] = Plot_Shadow_Manifold()
	N = 1000;
	t = linspace(0,10*3.14, N);
	X_temp = 0.75*sin(0.25*t) + 0.125*sin(20*t);

	plot(t,X_temp)
	xlim([0 10*3.14])
	uiwait

	Start = 100;
	Stop = 900;

	% the high oscillation frequency is 5 time indexes
	for tau = 1:2:15
		X = X_temp(Start:Stop);
		Y = X_temp(Start-tau:Stop-tau);
		Z = X_temp(Start-2*tau:Stop-2*tau);
	
		figure1 = figure('Name','Phase Space Example');
		axes1 = axes('Parent',figure1);
		hold(axes1,'on');
		plot3(X,Y,Z,'Parent',axes1)
		title('3D Shadow Manifold M_X')
		xlabel('X(t)') 
		ylabel('X(t-\tau)') 
		zlabel('X(t-2\tau)') 
		strang = strcat('tau = ',mat2str(tau));
		legend(strang)
		uiwait
	end
endfunction