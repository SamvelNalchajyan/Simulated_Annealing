#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

struct Condition
{
	double x;
	double y;
};

std::vector<std::vector<double>> A_matrix = { {-0.940, -0.536, -0.743}, {-0.502,  0.804,  0.769}, {-0.428, -0.789,  0.204} };
std::vector<std::vector<double>> B_matrix = { { 0.590,  0.160, -0.681}, { 0.387,  0.945, -0.195}, {-0.231,  0.152,  0.295} };
std::vector<std::vector<double>> C_matrix = { {-0.896, -0.613, -0.463}, { 0.038, -0.428, -0.714}, { 0.103,  0.741, -0.317} };
std::vector<std::vector<double>> D_matrix = { {-0.754, -0.558, -0.989}, {-0.702,  0.881,  0.397}, {-0.056,  0.085, -0.616} };

double horrific_function(std::vector<double> argument)
{
	if (argument.size() == 2)
	{
		double x = argument[0];
		double y = argument[1];
		double result = 0;
		double sum1 = 0, sum2 = 0;
		double A, B, C, D;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				A = A_matrix[i][j];
				B = B_matrix[i][j];
				C = C_matrix[i][j];
				D = D_matrix[i][j];
				sum1 += A * sin(i * M_PI * (x - 1 / 2)) * sin(j * M_PI * (y - 1 / 2)) + B * cos(i * M_PI * (x - 1 / 2)) * cos(j * M_PI * (y - 1 / 2));
				sum2 += C * sin(i * M_PI * (x - 1 / 2)) * sin(j * M_PI * (y - 1 / 2)) + D * cos(i * M_PI * (x - 1 / 2)) * cos(j * M_PI * (y - 1 / 2));
			}
		}
		return sqrt(sum1 * sum1 + sum2 * sum2);
	}
	else
	{
		throw - 1;
	}
}

double TestFunction(std::vector<double> argument)
{
	if (argument.size() == 2)
	{
		double x = argument[0];
		double y = argument[1];
		//return (1 + sin(10 * x) + cos(2 * x) + cos(2 * x + 2 * y) + cos(2 * y) + sin(20 * y) + y * y);
		return ((1 - x) * (1 - x) + 100 * (y - x * x) * (y - x * x));
	}
	else
	{
		throw - 1;
	}
}

//t_0 - начальная температура
//i - номер итерации
double Temperature(double t_0, int i)
{
	double k = 1; 
	return t_0 / (k * i);
}

//dE - разность целевой функции
//t - температура в текущей точке
bool Transition(double dE, double t)
{
	double P;
	if (dE < 0)
	{
		P = 1;
	}
	else
	{
		P = exp(-dE / t);
	}
	double rand_value = static_cast<double>(rand()) / RAND_MAX;
	if (rand_value <= P)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//curr - текущее состояние
std::vector<double> GenerateCandidate(std::vector<double> curr, double x_min, double x_max, double y_min, double y_max)
{
	double k_x = 0.000001;
	double k_y = 0.000001;
	std::vector<double> res;

	double x = curr[0] + (k_x * ((2 * rand()) - RAND_MAX) * (rand() % 11));
	double y = curr[1] + (k_y * ((2 * rand()) - RAND_MAX) * (rand() % 11));

	if (x >= x_min && x <= x_max)
	{
		res.push_back(x);
	}
	else
	{
		if (x < x_min)
		{
			x = x_min;
			res.push_back(x);
		}
		else
		{
			x = x_max;
			res.push_back(x);
		}
	}

	if (y >= y_min && y <= y_max)
	{
		res.push_back(y);
	}
	else
	{
		if (y < y_min)
		{
			y = y_min;
			res.push_back(y);
		}
		else
		{
			y = y_max;
			res.push_back(y);
		}
	}
	return res;
}

/*
Condition GenerateCandidate(Condition curr)
{

	double k_x = 0.00001;
	double k_y = 0.00001;
	Condition res;

	res.x = curr.x + k_x * ((2 * rand()) - RAND_MAX);
	res.y = curr.y + k_y * ((2 * rand()) - RAND_MAX);

	return res;
}
*/

//start - начальная точка
//initialTemperature - начальная температура
//endTemperature - конечная температура
std::vector<double> SimulatedAnnealing(std::vector<double> start, double initialTemperature, double endTemperature, double x_min, double x_max, double y_min, double y_max)
{
	double t = initialTemperature;
	std::vector<double> s_prev(start);
	std::vector<double> s_curr(start);
	std::vector<double> candidate;
	double dE;
	int i = 1;
	while (t > endTemperature)
	{
		candidate = GenerateCandidate(s_prev, x_min, x_max, y_min, y_max);
		//dE = TestFunction(candidate) - TestFunction(s_prev);
		dE = horrific_function(candidate) - horrific_function(s_prev);
		if (Transition(dE, t))
		{
			s_curr = candidate;
		}
		t = Temperature(initialTemperature, i);
		i++;
		s_prev = s_curr;
	}
	return s_curr;
}