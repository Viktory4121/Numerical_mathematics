// Вариант 2.
#define _USE_MATH_DEFINES
#include <iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<string>

using namespace std;

//---------------------ПО ВАРИАНТУ---------------------
const double Ax = 0;
const double Bx = 1;
const double At = 0;
const double Bt = 0.1;
const int A = 1;
//---------------------
//f(t)
double f_t(double x, double t, double coef = 1) {
	return coef * ((1 - t) * exp(-x));
}

//u(0,x) = 0
vector<double> u_0_x(vector<double> x) {
	vector<double> u(x.size(), 0);
	return u;
}

//u(t,0) = t
vector<double> u_t_0(vector<double> t) {
	vector<double> u;
	for (int i = t.size() - 1; i >= 1; i--) u.push_back(t[i]);
	return u;
}
double ut_0(double t) {
	return t;
}

//u(t,1) = t/exp
vector<double> u_t_1(vector<double> t) {
	vector<double> u;
	for (int i = t.size() - 1; i >= 1; i--) u.push_back(t[i] / M_E);
	return u;
}
double ut_1(double t) {
	return t / M_E;
}
//-----------------------------------------------------
//Генерация переменных x и t
void gen_var(double step_x, double step_t, vector<double> *x, vector<double> *t) {
	
	for (double i = Ax; i <= Bx; i += step_x) {
		(*x).push_back(i);
	}
	if ((round((*x)[(*x).size() - 1] * 10.0) / 10.0) != Bx) (*x).push_back(Bx);

	for (double i = At; i <= Bt; i += step_t) {
		(*t).push_back(i);
	}
	if ((round((*t)[(*t).size() - 1] * 1000.0) / 1000.0) != Bt) (*t).push_back(Bt);
}


// Метод прогонки, возвращающий строку
vector<double> progon(double h, double tau, int n, vector<vector<double>>* table, int index_T, double k = 1)
{
	vector<double> delta(n, 0), lambda(n,0), x(n,0),
		d (n, (-(tau * A * A) / (h * h))),
		c (n, (1 + (2 * tau * A * A) / (h * h))),
		b (n, (-(tau * A * A) / (h * h))),
		r (n, 0);

	for (int i = 0; i < n; i++) { //заполнение вектора r[j]
		if (i == 0) {
			r[i] = (*table)[index_T - 1][1] + (tau * k * f_t(h, tau * (index_T))) + (tau * A * A * ut_0(tau * (index_T)) / (h * h));
		} else if (i == n - 1) {
			r[i] = (*table)[index_T - 1][n] + (tau * k * f_t(n * h, (index_T) * tau)) + (tau * A * A * ut_1(tau * (index_T)) / (h * h));
		} else {
			r[i] = ((*table)[index_T - 1][i + 1] + (tau * k * f_t((i + 1) * h, (index_T) * tau)));
		}
	}
	
	//Вспомогательные величины:
	delta[0] = -d[0] / c[0];
	lambda[0] = r[0] / c[0];

	for (int i = 1; i < n; i++) {
		delta[i] = -d[i] / (c[i] + delta[i - 1] * b[i]);
		lambda[i] = (r[i] - b[i] * lambda[i - 1]) / (c[i] + delta[i - 1] * b[i]);
	}

	x[n - 1] = lambda[n - 1];

	for (int i = n - 2; i >= 0; i--) { //Расчет значений строки
		x[i] = delta[i] * x[i + 1] + lambda[i];
	}

	return x;
}

//Вывод в файл
void output_file(vector<vector<double>> u1, vector<vector<double>> u2, vector<double> x, vector<double> t, const char* name_file, const char* message) {
	FILE* file;
	fopen_s(&file, name_file, "w");

	fprintf_s(file, "Лабораторная работа №5. Вариант 2.\n");
	fprintf_s(file, "Численные методы решения ДУЧП параболического типа.");

	fprintf_s(file, "\nПолучили таблицу с коэффициентом при f(t) равном 1.\n");
	fprintf_s(file, message);
	int tt = t.size() - 1;
	for (int i = 0; i < t.size(); i++) {
		fprintf_s(file, "%7.4f|", t[tt]);							//вывод столбца t

		for (int j = 0; j < x.size(); j++) {
			fprintf_s(file, "%11.4f", u1[i][j]);						//вывод значений функции u	
		}
		fprintf_s(file, "\n");
		tt--;
	}
	fprintf_s(file, " t / x |");
	for (int i = 0; i < x.size(); i++) {
		fprintf_s(file, "%11.4f", x[i]);								//вывод строки x		
	}
	//---------------------------------------------------
	fprintf_s(file, "\n\nПолучили таблицу с коэффициентом при f(t) равном 1.1.\n");
	fprintf_s(file, message);
	int ttt = t.size() - 1;
	for (int i = 0; i < t.size(); i++) {
		fprintf_s(file, "%7.4f|", t[ttt]);								//вывод столбца t

		for (int j = 0; j < x.size(); j++) {
			fprintf_s(file, "%11.4f", u2[i][j]);						//вывод значений функции u	
		}
		fprintf_s(file, "\n");
		ttt--;
	}
	fprintf_s(file, " t / x |");
	for (int i = 0; i < x.size(); i++) {
		fprintf_s(file, "%11.4f", x[i]);								//вывод строки x		
	}

	fclose(file);
}

//Явный метод
void explicit_method(double h, double tau, bool logic) {
	vector<double> x, t;
	gen_var(h, tau, &x, &t);
	vector<vector<double>> u1(t.size(), vector<double>(x.size())), u2(t.size(), vector<double>(x.size()));
	vector<double> ut0 = u_t_0(t);
	vector<double> u0x = u_0_x(x);
	vector<double> ut1 = u_t_1(t);

	int ind_t = t.size() - 1;
	for (int i = 0; i < t.size() - 1; i++) {
		u1[i][0] = ut0[i];
		u1[i][x.size() - 1] = ut1[i];
		u2[i][0] = ut0[i];
		u2[i][x.size() - 1] = ut1[i];
	}
	for (int i = 0; i < x.size(); i++) {
		u1[t.size() - 1][i] = u0x[i];
		u2[t.size() - 1][i] = u0x[i];
	}


	for (int i = t.size() - 2; i >= 0; i--) {
		for (int j = 1; j <= x.size() - 2; j++) {
			u1[i][j] = A * A * tau *((u1[i + 1][j + 1] - 2 * u1[i + 1][j] + u1[i + 1][j - 1]) / (h * h)) + tau * f_t(x[j - 1], t[i + 1]) + u1[i + 1][j];
			u2[i][j] = A * A * tau * ((u2[i + 1][j + 1] - 2 * u2[i + 1][j] + u2[i + 1][j - 1]) / (h * h)) + tau * f_t(x[j - 1], t[i + 1], 1.1) + u2[i + 1][j];
		}
	}

	if (logic) {
		output_file(u1, u2, x, t, "Явный метод с выполненным условием Куранта.txt", "Явный метод с выполненным условием Куранта выглядит так:\n");
	} else {
		output_file(u1, u2, x, t, "Явный метод с невыполненным условием Куранта.txt", "Явный метод с невыполненным условием Куранта выглядит так:\n");
	}

	//файл под названием "Явный метод с выполненным условием Куранта.txt" будет содержать 2 таблицы для coef=1 и coef=1,1 
	//файл под названием "Явный метод с не выполненным условием Куранта.txt" будет содержать 2 таблицы для coef=1 и coef=1,1
}

//Неявный метод
void implicit_method(double h, double tau, bool logic) {
	vector<double> x, t;
	gen_var(h, tau, &x, &t);
	vector<vector<double>> u1(t.size(), vector<double>(x.size())), u2(t.size(), vector<double>(x.size()));
	vector<double> ut0 = u_t_0(t);
	vector<double> u0x = u_0_x(x);
	vector<double> ut1 = u_t_1(t);

	for (int i = 0; i < t.size() - 1; i++) {
		u1[i][0] = ut0[i];
		u1[i][x.size() - 1] = ut1[i];
		u2[i][0] = ut0[i];
		u2[i][x.size() - 1] = ut1[i];
	}
	for (int i = 0; i < x.size(); i++) {
		u1[t.size() - 1][i] = u0x[i];
		u2[t.size() - 1][i] = u0x[i];
	}

	vector<double> s1, s2;
	for (int i = 1; i <= t.size() - 1; i++) {
		s1 = progon(h, tau, x.size() - 2, &u1, i);
		s2 = progon(h, tau, x.size() - 2, &u2, i, 1.1);
		
		for (int j = 0; j < s1.size(); j++) {
			u1[t.size() - i - 1][j + 1] = s1[j];
			u2[t.size() - i - 1][j + 1] = s2[j];
		}
	}
	

	if (logic) {
		output_file(u1, u2, x, t, "Неявный метод с выполненным условием Куранта.txt", "Неявный метод с выполненным условием Куранта выглядит так:\n");
	}
	else {
		output_file(u1, u2, x, t, "Неявный метод с невыполненным условием Куранта.txt", "Неявный метод с невыполненным условием Куранта выглядит так:\n");
	}
	//файл под названием "Неявный метод с выполненным условием Куранта.txt" будет содержать 2 таблицы для coef=1 и coef=1,1
	//файл под названием "Неявный метод с не выполненным условием Куранта.txt" будет содержать 2 таблицы для coef=1 и coef=1,1
}


void calculating(double h, double tau, bool logic) {
	explicit_method(h, tau, logic);
	implicit_method(h, tau, logic);
}


int main() {
	setlocale(LC_ALL, "Russian");

	double tau;
	double h;
	bool logic = false;

	cout << "Введите шаг h для x: ";
	cin >> h;
	cout << "Введите шаг tau для t: ";
	cin >> tau;

	if ((A * A * (tau / (h * h))) <= 0.5) {
		cout << "\nВыполнено условие Куранта.\n";
		logic = true;
	}
	else {
		cout << "\nУсловие Куранта не выполнено.\n";
	}

	calculating(h, tau, logic);

	cout << "\nПроизведёны выводы в файлы с результатом выполнения программы.";

	return 0;
}