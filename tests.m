% Thomas Merkh, tmerkh@ucla.edu, Last Update: July 25th, 2017
% This test script runs the Conv_Cross_Map program to determine causation between two time series
% Attempts to determine whether or not X(t) causes Y(t) and/or Y(t) causes X(t)

% Things to try - randn noise, and effect of shadow dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test Data Sizes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STARTTIME = 50;							% Indicates where to begin constructing time series library. Require STARTTIME > 2*E or knn search will throw an indexing error.
TIME = 1000;							% Typically indicates the final time of the time series data
number_time_points = 10;				% Size of time series library - points uniformly spaced over [STARTTIME, TIME]
Intv = floor(TIME/number_time_points);	% Spacing between time series points in library
tau = 2;								% Amount of time shift in reconstructing the phase space - Must be Integer!!
plot_phase_space = true;				% Plot 3D shadow manifold of phase space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Test Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test# should be a 'TIME by 2' matrix containing two time series that may or may not be causally defined.
time_vec = linspace(0,2*3.14,TIME);
k = floor(TIME*0.8); 			
noise_strength = 0.8;
noise = zeros(1,TIME);
noise(randi([1 TIME], 1, k)) = randn(1,k);
noise2 = zeros(1,TIME);
noise2(randi([1 TIME], 1, k)) = randn(1,k);
noise3 = zeros(1,TIME);
noise3(randi([1 TIME], 1, k)) = randn(1,k);
Zt_noise = 2*time_vec.*cos(time_vec) + noise*noise_strength;  % Used to muddle the correlations
Zt_noise(1:floor(TIME/3)) = 0;								  % Very dramatic comparison - Turns on for t > T/3

% Test1
Yt = sin(time_vec);  						% Second Time Series
Xt = Yt + 3 + Zt_noise;						% First Time Series
test1 = [Xt', Yt'];

% Test2
Yt = 3*exp(-1*time_vec.*time_vec);  		
Xt = 2*Yt.*cos(3*time_vec);					
test2 = [Xt', Yt'];

% Test3
Yt = exp(-1*time_vec).*sin(5*time_vec);    
Xt = 4*Yt.*Yt.*sin(time_vec) - 1;		
test3 = [Xt', Yt'];

% Test 4
Yt = 10*exp(-1*time_vec).*sin(time_vec) + 3*exp(-1*(time_vec-5.14).*(time_vec-5.14)).*cos(10*time_vec) + noise2*noise_strength;
Xt = -1*Yt.*Yt + Yt - Zt_noise - noise3*noise_strength;													 
test4 = [Xt', Yt'];

% Test 5 - ICS set up on the brink of instability
Xt = zeros(1,TIME); Xt(1) = 0.585;
Yt = zeros(1,TIME); Yt(1) = 1.10;
dt = 2*3.14/TIME;
for i = 1:TIME-1
	Xt(i+1) = Xt(i) + dt*(Xt(i)*Yt(i) - Yt(i)*Yt(i));
	Yt(i+1) = Yt(i) - dt*(Xt(i)*Yt(i) + Xt(i)*Xt(i));
end
Xt = Xt + 0.1*noise2;
Yt = Yt + 0.1*noise;
test5 = [Xt', Yt'];

% Test 6 - Smooth multiscale time series
time_vec = linspace(0,10*3.1415, TIME);
Xt = sin(0.25*time_vec);
Yt = sin(21.5*time_vec);
Zt = Xt + Yt;
test6 = [Xt', Yt', Zt'];

% Show Test Case
figure0_1 = figure('Name','Test Case');
axes0_1 = axes('Parent',figure0_1);
hold(axes0_1,'on');
plot(time_vec,test6(:,2),time_vec,test6(:,3),'Parent',axes0_1);
title('Time Series')
xlabel('Time') 
ylabel('Measurements') 
xlim([0 max(time_vec)])
legend('X','Y')

if(plot_phase_space)
	Z1 = Zt(STARTTIME:TIME);
	Z2 = Zt(STARTTIME-tau:TIME-tau);
	Z3 = Zt(STARTTIME-2*tau:TIME-2*tau);
	
	figure0_2 = figure('Name','3D Shadow Manifold M_Y of Test Case');
	axes0_2 = axes('Parent',figure0_2);
	hold(axes0_2,'on');
	plot3(Z1,Z2,Z3,'Parent',axes0_2,'k','LineWidth',2)
	title('3D Shadow Manifold M_Y')
	xlabel('Y(t)') 
	ylabel('Y(t-\tau)') 
	zlabel('Y(t-2\tau)') 
	strang = strcat('tau = ',mat2str(tau));
	legend(strang)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin CCM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = test6(:,2);
Y = test6(:,3);
E = 3;              			% Embedding dimension for reconstruction - Usually a low number (3, ..., 6) is desirable. 
L = STARTTIME:Intv:TIME;        % Points in time series in which one would like to compute causation at.  Typically, uniformly spread out is ideal.

Corr_Coefficients = zeros(2,length(L));
for i = 1:length(L)
    [Corr_Coefficients(:,i), Predicted_Y, Predicted_X, original_Y, original_X] = Conv_Cross_Map(X(1:L(i)),Y(1:L(i)),tau,E); 
end

% Plot original data and estimated data using CCM
% Y is such a scattered plot because the cross mapping of Y using M_x fails - x doesn't predict y well!!
% X is a nice linear looking plot because one can accurately predict C using M_y!

figure1 = figure('Name','The Data');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(original_X,Predicted_X,'o',original_Y,Predicted_Y,'*','Parent',axes1)
title('Estimated vs. original data')
xlabel('Original data') 
ylabel('Estimated data') 
legend('X','Y')

% Plot Degree to which X predicts Y and Y predicts X
% X|M_Y strong correlation means X depends on Y strongly!
figure2 = figure('Name','CCM Correlation Strength');
axes2 = axes('Parent',figure2);
hold(axes2,'on');
plot(L,Corr_Coefficients,'Parent',axes2, 'LineWidth', 2)
leg = legend('X(t)|M_Y','Y(t)|M_X');  
set(leg, 'Interpreter', 'latex');
xlabel('Time','Interpreter', 'latex')
ylabel('Correlation Coefficient','Interpreter', 'latex')
title('X(t) dependence on Y(t), and vise versa')