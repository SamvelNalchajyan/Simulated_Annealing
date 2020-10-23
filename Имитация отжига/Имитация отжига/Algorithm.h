#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

struct Condition
{
	double x;
	double y;
};

double TestFunction(Condition s)
{
	//return (1 + sin(10 * s.x) + cos(2 * s.x) + cos(2 * s.x + 2 * s.y) + cos(2 * s.y) + sin(20 * s.y) + s.y * s.y);
	return ((1 - s.x) * (1 - s.x) + 100 * (s.y - s.x * s.x) * (s.y - s.x * s.x));
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

Condition GenerateCandidate(Condition curr)
{
	double k_x = 0.000001;
	double k_y = 0.000001;
	Condition res;

	res.x = curr.x + (k_x * ((2 * rand()) - RAND_MAX) * (rand() % 11));
	res.y = curr.y + (k_y * ((2 * rand()) - RAND_MAX) * (rand() % 11));

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
Condition SimulatedAnnealing(Condition start, double initialTemperature, double endTemperature)
{
	double t = initialTemperature;
	Condition s_prev = start;
	Condition s_curr = start;
	Condition candidate;
	double dE;
	int i = 1;
	while (t > endTemperature)
	{
		candidate = GenerateCandidate(s_prev);
		dE = TestFunction(candidate) - TestFunction(s_prev);
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