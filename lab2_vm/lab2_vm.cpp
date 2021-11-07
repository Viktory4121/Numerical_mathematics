//Вариант 2. f(x)=pow(tan(x), 2)+pow(1/tan(x),2); [pi/6, pi/3]
//Численное интегрирование (формулы трапеций, Симпсона, Ньютона)

#define _USE_MATH_DEFINES
#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<ctime>

//----------ПО---ВАРИАНТУ----------
const double A = M_PI / 6;
const double B = M_PI / 3;
//---------------------------------

using namespace std;

//----------ПО---ВАРИАНТУ----------
//вызов функции (возвращает 1 элемент)
double function(double x) {
    return (pow(tan(x), 2) + pow((1 / tan(x)), 2));
}
//---------------------------------

//дублируется вывод в файл
void Output_file(double epsilon, double trapezoid, double second_T, double simpson, double second_S, double newton, double second_N,
				 double n_T, double n_S, double n_N) {

	ofstream file;
	file.open("output.txt");

	file << "Введённое значение epsilon = " << epsilon << endl;
	file << "Вычисление приближённого значения определённого интеграла" << endl;
	file << "Вариант 2. f(x)=tan(x)^2 + (1/tan(x))^2; [pi/6, pi/3]" << endl << endl;
	file << fixed << setprecision(6) << "По методу трапеций: " << trapezoid << "; n = " << int(n_T) << endl;
	file << fixed << setprecision(6) << "По методу Симпсона: " << simpson << "; n = " << int(n_S) << endl;
	file << fixed << setprecision(6) << "По методу Ньютона: " << newton << "; n = " << int(n_N) << endl << endl;

	file << "Истинное значение интеграла: " << 1.26220352556191 << endl << endl;
	
	file << fixed << setprecision(3) << "Время работы метода трапеций: в милисекундах - " << second_T << endl;
	file << fixed << setprecision(3) << "Время работы метода Симпсона: в милисекундах - " << second_S << endl;
	file << fixed << setprecision(3) << "Время работы метода Ньютона: в милисекундах - " << second_N << endl;

	file.close();
}

//консольный вывод
void Output_console(double trapezoid, double second_T, double simpson, double second_S, double newton, double second_N,
					double n_T, double n_S, double n_N) {

	cout << "Вычисление приближённого значения определённого интеграла" << endl;
	cout << "Вариант 2. f(x)=tan(x)^2 + (1/tan(x))^2; [pi/6, pi/3]" << endl << endl;
	cout << fixed << setprecision(5) << "По методу трапеций: " << trapezoid << "; n = " << int(n_T) << endl;
	cout << fixed << setprecision(5) << "По методу Симпсона: " << simpson << "; n = " << int(n_S) << endl;
	cout << fixed << setprecision(5) << "По методу Ньютона: " << newton << "; n = " << int(n_N) << endl << endl;

	cout << "Истинное значение интеграла: " << 1.26220352556191 << endl << endl;

	cout << fixed << setprecision(3) << "Время работы метода трапеций: в милисекундах - " << second_T << endl;
	cout << fixed << setprecision(3) << "Время работы метода Симпсона: в милисекундах - " << second_S << endl;
	cout << fixed << setprecision(3) << "Время работы метода Ньютона: в милисекундах - " << second_N << endl;
}

//-----------------------------------------------------------------------------
//расчет точности
double accuracy(int p, double int1, double int2) {
	return (abs((int1 - int2) / (pow(2, p) - 1)));
}

//-----------------------------------------------------------------------------
//метод трапеций для вычисления интеграла
double integral_trapezoid(double h) {
	vector<double> x;
	for (double i = A; i <= B; i += h)
		x.push_back(i);

	double integ = 0;
	for (int i = 1; i < x.size(); i++)
		integ += (x[i] - x[i - 1]) * (0.5 * function(x[i - 1]) + 0.5 * function(x[i]));

	return integ;
}

double Method_trapezoid(double eps, double *n) {
	double h = B - A, int1, int2, acc;
	
	int1 = integral_trapezoid(h);
	while (true) {
		int2 = integral_trapezoid(h/2);
		acc = accuracy(2, int1, int2);

		if (eps > acc) {
			*n = (B - A) / (h / 2); 
			return int2;
		}

		int1 = int2;
		h /= 2;
	}
}

//-----------------------------------------------------------------------------
//метод Симпсона для выч-я инт-а
double integral_simpson(double h) {
	vector<double> x;
	for (double i = A; i <= B; i += h)
		x.push_back(i);

	double integ = 0; //сюда будет суммироваться значение интеграла

	for (int i = 1; i < x.size(); i++)
		integ += (x[i] - x[i - 1]) * ((1.0 / 6.0) * function(x[i]) + (2.0 / 3.0) * function((x[i] + x[i - 1]) / 2.0) + (1.0 / 6.0) * function(x[i - 1]));

	return integ;
}

double Method_simpson(double eps, double *n) {
	double h = B - A, int1, int2, acc;

	int1 = integral_simpson(h);
	while (true) {
		int2 = integral_simpson(h / 2);

		acc = accuracy(4, int1, int2);
		//cout << int1 << " " << int2 << " " << acc << endl;
		if (eps > acc) {
			*n = (B - A) / (h / 2);
			return int2;
		}

		int1 = int2;
		h /= 2;
	}
}

//-----------------------------------------------------------------------------
//метод Ньютона для выч-я инт-ла

double integral_newton(double h) {
	vector<double> x;
	for (double i = A; i <= B; i += h)
		x.push_back(i);

	double integ = 0;
	for (int i = 1; i < x.size(); i++) {
		integ += (x[i] - x[i - 1]) * (0.125 * function(x[i - 1]) + 0.375 * function((2 * x[i - 1] + x[i]) / 3) 
			  + 0.375 * function((x[i - 1] + 2 * x[i]) / 3) + 0.125 * function(x[i]));
	}

	return integ;
}

double Method_newton(double eps, double *n) {
	double h = B - A, int1, int2, acc;

	int1 = integral_newton(h);
	while (true) {
		int2 = integral_newton(h / 2);
		acc = accuracy(4, int1, int2);

		if (eps > acc) {
			*n = (B - A) / (h / 2);
			return int2;
		}

		int1 = int2;
		h /= 2;
	}

}

//-----------------------------------------------------------------------------
//вызов функций для расчёта интеграла разными методами и замер времени работы методов
void calculating_the_integral(double eps) {
	double n_T, n_S, n_N;

	double start_time_T = clock();
	double trap = Method_trapezoid(eps, &n_T);
	double end_time_T = clock();
	double second_T = end_time_T - start_time_T;

	double start_time_S = clock();
	double simp = Method_simpson(eps, &n_S);
	double end_time_S = clock();
	double second_S = end_time_S - start_time_S;

	double start_time_N = clock();
	double newt = Method_newton(eps, &n_N);
	double end_time_N = clock();
	double second_N = end_time_N - start_time_N;

	Output_console(trap, second_T, simp, second_S, newt, second_N, n_T, n_S, n_N);
	Output_file(eps, trap, second_T, simp, second_S, newt, second_N, n_T, n_S, n_N);
}

int main() {
    setlocale(LC_ALL, "Russian");

	double eps_0;

	cout << "Введите значение epsilon = ";
	cin >> eps_0;

	calculating_the_integral(eps_0);

	cout << endl << "Запись в файл продублирована." << endl;

    return 0;
}