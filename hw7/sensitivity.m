%% Parameters
NUM_POINTS = 200;
SAVE = true;
LOWER_LIM = -3;
UPPER_LIM = 5;

%% Compile Fortran Code
system('make clean');
system('make');

%% Setup
u_truth = readmatrix("output/truth.dat");

% Initialize vectors
model_data_misfits = logspace(LOWER_LIM, UPPER_LIM, NUM_POINTS);
rms = zeros(size(model_data_misfits));

%% Perform data inversion for each model-data misfit value
for i=1:NUM_POINTS
	disp([num2str(i) ' of ' num2str(NUM_POINTS)]);
	covariance_matrices(model_data_misfits(i));
	system('./invert_data');
	
	u_fit = readmatrix("output/u_fit.dat")';
	
	error = u_fit - u_truth;
	
	rms(i) = sqrt(mean(error.^2));
end

%% Reset outputs for other analysis scripts
covariance_matrices;
system('./invert_data');

%% Get RMS from original value
u_fit = readmatrix("output/u_fit.dat")';
error = u_fit - u_truth;
rms_og = sqrt(mean(error.^2));

%% Plotting
close all;

semilogx(model_data_misfits, rms);
hold on;
plot(0.05, rms_og, 'x');
hold off;
title("Inverse Bias vs Model-Data Misfit");
xlabel("Model-Data Misfit");
ylabel("RMS of Inverse Bias");
xlim([10^LOWER_LIM, 10^UPPER_LIM]);
legend("RMS Values", "Standard RMS Value", "Location", "northwest");

if SAVE
	saveas(gcf, 'figures/inverse_bias_vs_model-data_misfit.png');
end
