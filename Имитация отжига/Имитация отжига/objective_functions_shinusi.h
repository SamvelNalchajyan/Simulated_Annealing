#pragma once
#include <vector>
#define PI 3.14159265358979323846

double Integral(double function(double, std::vector<double>), std::vector<double> argument) {
	int n = 20;
	const double width = 1.0 / (double)n;

	double simpson_integral = 0;
	for (int step = 0; step < n; step++) {
		const double t1 = step * width;
		const double t2 = (step + 1) * width;

		simpson_integral += (t2 - t1) / 6.0 * (function(t1, argument) + 4.0 * function(0.5 * (t1 + t2), argument) + function(t2, argument));
	}

	return simpson_integral;
}

double plancton_strategy(double t, std::vector<double> argument) {
	double result = argument[0];
	for (int i = 1; i < argument.size(); i++) {
		if (i % 2) {
			result += argument[i] * cos(2 * PI * ceil(i / 2.0) * t);
		}
		else {
			result += argument[i] * sin(2 * PI * ceil(i / 2.0) * t);
		}
	}
	return result;
}
double kinetic_energy_function(double t, std::vector<double> argument) {
	double result = 0;
	for (int i = 1; i < argument.size(); i++) {
		if (i % 2) {
			result -= argument[i] * 2 * PI * ceil(i / 2.0) * sin(2 * PI * ceil(i / 2.0) * t);
		}
		else {
			result += argument[i] * 2 * PI * ceil(i / 2.0) * cos(2 * PI * ceil(i / 2.0) * t);
		}
	}
	return result * result;
}

double food_function(double t, std::vector<double> argument) {
	const double sigma_1 = 0.25;
	const double C = 140;
	const double C1 = 40;

	double x = plancton_strategy(t, argument);
	double result;
	if ((x > -C) && (x < 0))
		result = tanh(sigma_1 * (x + C1)) + 1;
	else
		result = 0;
	return result;
}

double predator_depth_function(double t, std::vector<double> argument) {
	const double sigma_2 = 0.003;
	const double C = 140;
	const double C1 = 40;

	double x = plancton_strategy(t, argument);
	double result;
	if ((x > -C) && (x < 0))
		result = tanh(sigma_2 * (x + C1)) + 1;
	else
		result = 0;
	return result;
}

double predator_time_function(double t) {
	const double epsilon = 0.013;
	double res = 0;
	if (t > 0 && t < 1)
		res = cos(2 * PI * t) - epsilon * cos(6 * PI * t) + 1;
	else
		res = 0;
	return res;
}

double predator_function(double t, std::vector<double> argument) {
	return predator_depth_function(t, argument) * predator_time_function(t);
}

double mortality_function(double t, std::vector<double> argument) {
	const double C0 = 60;
	const double ksi = 5e-19;

	return ksi * cosh(plancton_strategy(t, argument) + C0);
}

double rosenbrock_function(std::vector<double> argument) {
	if (argument.size() == 2) {
		double x = argument[0];
		double y = argument[1];
		return pow((1 - x), 2) + 100 * pow((y - x * x), 2);
	}
	else
		throw - 1;
}
double icicle_function(std::vector<double> argument) {
	if (argument.size() == 2) {
		double x = argument[0];
		double y = argument[1];
		return (1 + sin(10 * x) + cos(2 * x) + cos(2 * x + 2 * y) + cos(2 * y) + sin(20 * y) + y * y);
	}
	else
		throw - 1;

}
double horrific_function(std::vector<double> argument) {
	if (argument.size() == 2) {
		double x = argument[0];
		double y = argument[1];
		double result = 0;
		double sum1 = 0, sum2 = 0;
		double A, B, C, D;

		std::vector<std::vector<double>> A_matrix =
		{
													  {-0.940, -0.536, -0.743},
													  {-0.502,  0.804,  0.769},
													  {-0.428, -0.789,  0.204}
		};
		std::vector<std::vector<double>> B_matrix =
		{
													  { 0.590,  0.160, -0.681},
													  { 0.387,  0.945, -0.195},
													  {-0.231,  0.152,  0.295}
		};
		std::vector<std::vector<double>> C_matrix =
		{
													  {-0.896, -0.613, -0.463},
													  { 0.038, -0.428, -0.714},
													  { 0.103,  0.741, -0.317}
		};
		std::vector<std::vector<double>> D_matrix =
		{
													  {-0.754, -0.558, -0.989},
													  {-0.702,  0.881,  0.397},
													  {-0.056,  0.085, -0.616}
		};

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				A = A_matrix[i][j];
				B = B_matrix[i][j];
				C = C_matrix[i][j];
				D = D_matrix[i][j];
				sum1 += A * sin(i * PI * (x - 1 / 2)) * sin(j * PI * (y - 1 / 2)) + B * cos(i * PI * (x - 1 / 2)) * cos(j * PI * (y - 1 / 2));
				sum2 += C * sin(i * PI * (x - 1 / 2)) * sin(j * PI * (y - 1 / 2)) + D * cos(i * PI * (x - 1 / 2)) * cos(j * PI * (y - 1 / 2));
			}
		}
		return sqrt(sum1 * sum1 + sum2 * sum2);
	}
	else
		throw - 1;
}
double plancton_function(std::vector<double> argument) {
	const double C = 140;
	const double alpha = 2;
	const double beta = 2.5e-5;
	const double gamma = 333;
	const double delta = 0.01;

	double result = 0;
	double I1 = alpha * Integral(food_function, argument);
	result += I1;

	double I2 = beta * Integral(kinetic_energy_function, argument);
	result -= I2;

	double I3 = gamma * Integral(predator_function, argument);
	result -= I3;
	
	double I4 = delta * Integral(mortality_function, argument);
	result -= I4;

	//std::cout << "(" << I1 << ", " << I2 << ", " << I3 << ", " << I4 << ")" << std::endl;;

	return result * (-1);
}
