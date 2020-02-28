% Thomas Merkh, tmerkh@ucla.edu, Updated Last: July 25th, 2017
% This function does the convergent cross mapping between two time series X and Y.
% E is the dimension of the shadow manifolds (topographically the same as the phase space manifold) i.e embedding dimension
% tau is the time shift used to construct the shadow manifolds, M_X, M_Y, of dimension E (Taken's Theorem)

% Output - OrigX, OrigY are just the original time series
% Predicted_X and Predicted_Y are the predicted (estimated) data.  
% If Y is accuractely predictable using M_X, then Predicted_X vs. original_X should be a correlated plot
% If X is not accurately predicted by M_Y, then Predicted_Y vs. original_Y should be an uncorrelated plot because CCM doesn't predict it well.
% Predicted_Corr is the correlation coefficients as a function of time

function [Predicted_Corr, Predicted_Y, Predicted_X, original_Y, original_X] = Conv_Cross_Map(X, Y, tau, E)

L = length(X);
T = 1+(E-1)*tau;  
% Largest time shift (increment up by tau, E times to get to this time). 
% When tau = 1, this is the number of time shifts - Equals the number of dimensions
N = L-T+1;			% Number of points on shadow manifold we will do the nearest neighbor search at
% Matrix containing the value of the time series X or Y with no shift, 1 shift, ...., 
% up to the number of shifts corresponding to dim E, i.e dimension of reconsturcted M_X or M_Y
% It does this for N points, which are our measurements on the manifold 
Xm = zeros(N,E);			
Ym = zeros(N,E);
Predicted_N = E+1;

%% Reconstructions of the original systems - "Library of points from each time series"
for t = 1:N
    Xm(t,:) = X((T+t-1):-tau:(T+t-1-(E-1)*tau));
    Ym(t,:) = Y((T+t-1):-tau:(T+t-1-(E-1)*tau));
end

Predicted_X = zeros(N,1);
Predicted_Y = zeros(N,1);

original_Y = Y(T:end);
original_X = X(T:end);

%% Neighborhood search - find E+2 nearest neighbors to each point of interest (in Xm)
% 'parfor' means each step can be done in parallel
parfor j = 1:N

	% Returns the neighborIds neighborDistances to the nearest E+2 points
	[n1,d1] = knnsearch(Xm,Xm(j,:),E+2);
	[n2,d2] = knnsearch(Ym,Ym(j,:),E+2);

	%% CMM - Algorithm Cited in Reference Paper Supplementary Materials
	% Causality studied in reconstructed state space,
	% Examples of uni-directionally connected chaotic systems, by
	% Anna Krakovská, Jozef Jakubík, Hana Budáčová, Mária Holecyová
	Temp_Y = original_Y(n1(2:end));
	Temp_X = original_X(n2(2:end));
	Predict_Temp_Y = Temp_Y(1:Predicted_N);
	Predict_Temp_X = Temp_X(1:Predicted_N);
	Predicted_d1 = d1(:,2:Predicted_N+1);
	Predicted_d2 = d2(:,2:Predicted_N+1);
	u1 = exp(-Predicted_d1./(Predicted_d1(:,1)*ones(1,Predicted_N)));
	u2 = exp(-Predicted_d2./(Predicted_d2(:,1)*ones(1,Predicted_N)));
	w1 = u1./(sum(u1,2)*ones(1,Predicted_N));
	w2 = u2./(sum(u2,2)*ones(1,Predicted_N));
	Predicted_Y(j) = w1*Predict_Temp_Y;
	Predicted_X(j) = w2*Predict_Temp_X;
end

Predicted_Corr1 = corrcoef(original_Y,Predicted_Y);
Predicted_Corr(2,1) = Predicted_Corr1(1,2);
Predicted_Corr2 = corrcoef(original_X,Predicted_X);
Predicted_Corr(1,1) = Predicted_Corr2(1,2);
end