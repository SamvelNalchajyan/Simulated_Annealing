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
	std::cout << "������� ��������� �����" << std::endl;
	Condition start;
	std::cin >> start.x >> start.y;
	std::cout << std::endl << "������� ��������� � �������� �����������" << std::endl;
	double t_max, t_min;
	std::cin >> t_max >> t_min;
	Condition res;
	for (int i = 0; i < 15; i++)
	{
		res = SimulatedAnnealing(start, t_max, t_min);
		std::cout << std::endl << "����� - (" << res.x << ", " << res.y << ")" << std::endl;
	}
	/*-------------------------------------------------*/
	system("pause");
	return 0;
}