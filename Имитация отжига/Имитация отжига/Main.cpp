#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <clocale>
#include "Algorithm.h"

int main()
{
	setlocale(LC_CTYPE, "Russian");
	/*-------------------------------------------------*/
	std::cout << "¬ведите начальную точку" << std::endl;
	std::vector<double> start;
	double x0, y0;
	std::cin >> x0 >> y0;
	start.push_back(x0);
	start.push_back(y0);
	//std::cout << std::endl << start.size() << std::endl;
	std::cout << std::endl << "¬ведите границы аргументов" << std::endl;
	double x_min, x_max, y_min, y_max;
	std::cin >> x_min >> x_max >> y_min >> y_max;
	std::cout << std::endl << "¬ведите начальную и конечную температуру" << std::endl;
	double t_max, t_min;
	std::cin >> t_max >> t_min;
	std::vector<double> res;
	for (int i = 0; i < 10; i++)
	{
		//start[0] = FRand(x_min, x_max);
		//start[1] = FRand(y_min, y_max);
		res = SimulatedAnnealing(start, t_max, t_min, x_min, x_max, y_min, y_max);
		std::cout << std::endl << "ќтвет - (" << res[0] << ", " << res[1] << ")" << std::endl;
	}
	std::vector<double> v1 = { -1.41, -0.39 };
	std::vector<double> v2 = { 1.73, -0.39 };
	std::cout << TestFunction(v1) << std::endl << TestFunction(v2) << std::endl;
	/*-------------------------------------------------*/
	system("pause");
	return 0;
}