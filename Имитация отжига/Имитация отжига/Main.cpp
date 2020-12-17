#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <clocale>
#include <ctime>
#include "Algorithm.h"

int main()
{
	setlocale(LC_CTYPE, "Russian");
	/*-------------------------------------------------*/
	std::vector<double> start;
	double x0 = 0, y0 = 0, z0 = 0, k0 = 0;
	start.push_back(x0);
	start.push_back(y0);
	//start.push_back(z0);
	//start.push_back(k0);
	
	//std::cout << "¬ведите начальную точку" << std::endl;
	//std::cin >> x0 >> y0;
	//start.push_back(x0);
	//start.push_back(y0);
	
	//std::cout << std::endl << start.size() << std::endl;
	std::cout << "¬ведите границы аргументов" << std::endl;
	double x_min, x_max, y_min, y_max;
	double z_min, z_max, k_min, k_max;
	std::cin >> x_min >> x_max >> y_min >> y_max;
	//std::cin >> z_min >> z_max >> k_min >> k_max;
	std::cout << std::endl << "¬ведите начальную и конечную температуру" << std::endl;
	double t_max, t_min;
	std::cin >> t_max >> t_min;
	std::vector<double> res;
	clock_t start_time, end_time;
	for (int i = 0; i < 10; i++)
	{
		start[0] = FRand(x_min, x_max);
		start[1] = FRand(y_min, y_max);
		//start[2] = FRand(y_min, y_max);
		//start[3] = FRand(y_min, y_max);
		start_time = clock();
		res = SimulatedAnnealing(start, t_max, t_min, x_min, x_max, y_min, y_max);
		end_time = clock();
		double duration = (double)(end_time - start_time) / CLOCKS_PER_SEC;
		std::cout << std::endl << "ќтвет - (" << res[0] << ", " << res[1] << ")" << std::endl << "«начение функции:" << plancton_function(res) << std::endl << "¬рем€ работы: " << duration << std::endl;
	}
	std::vector<double> v1 = { -85, -45 };
	//std::vector<double> v2 = { 1.73, -0.39 };
	//std::cout << plancton_function(v1) << std::endl;
	/*-------------------------------------------------*/
	system("pause");
	return 0;
}