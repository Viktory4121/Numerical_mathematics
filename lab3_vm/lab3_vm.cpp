//Вариант 2. Лабораторная работа №3.
//y' = f(x, y) = 2 * (2 - x) * y + 0.01 * exp(-x * x)
//x = [1, 6]
//y(А) = 10


#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<string>
#include <conio.h>

using namespace std;

//-------ПО ВАРИАНТУ--------------------
const double Y_A = 10.0; //значение у(A)
const double A = 1.0;	 //левая граница х
const double B = 6.0;	 //правая граница х

double function(double x, double y) {
	return (2 * (2 - x) * y + 0.01 * exp(-x * x));
}
//--------------------------------------
//расчет точности
double accuracy(int p, vector<double> y1, vector<double> y2) {
	vector<double> vv;
	double max_ = -1.0;
	int size_y;

	if (2 * y1.size() > y2.size()) { //проверка на выход за границы вектора
		size_y = y1.size() - 1;
	} else {
		size_y = y1.size();
	}

	for (int i = 0; i < size_y; i++) {
		vv.push_back(abs((y1[i] - y2[2 * i]) / (pow(2, p) - 1))); //ТУТ ВЫВАЛИВАЕТСЯ ИСКЛЮЧЕНИЕ-------------------------------------------------
		if (max_ < vv[i]) max_ = vv[i]; //поиск максимального значения в векторе
	}

	vv.clear();

	return max_;
}
//генерация х
vector<double> gen_x(double h) {
	vector<double> x;
	double i = A;

	while (i <= B) {
		x.push_back(i);
		i += h;
	}
	if (x[x.size() - 1] != B) x.push_back(B);

	return x;
}
//--------------------------------------
//МЕТОДЫ РЕШЕНИЯ:

vector<double> Method_Euler(double n, double eps, double *h_0) { //Метод Эйлера
	vector<double> y1, y2, x1, x2;
	double h0 = (B - A) / n;

	//вычисление для h
	x1 = gen_x(h0);
	y1.push_back(Y_A);
	for (int i = 1; i < x1.size(); i++) y1.push_back(y1[i - 1] + h0 * function(x1[i - 1], y1[i - 1]));
	//----------------

	while(true) {

		n *= 2;
		h0 = (B - A) / n;

		//вычисление для h/2
		x2 = gen_x(h0);
		y2.push_back(Y_A);
		for (int i = 1; i < x2.size(); i++) y2.push_back(y2[i - 1] + h0 * function(x2[i - 1], y2[i - 1]));
		//------------------

		double acc = accuracy(1, y1, y2); //порядок точности - 1

		if (eps > acc) { //проверка на точность
			*h_0 = h0; //чтобы запомнить точность
			return y2;
		}

		y1.clear();
		x1.clear();

		x1 = x2;
		y1 = y2;

		x2.clear();
		y2.clear();
	}
}

vector<double> Method_Hoina(double n, double eps, double* h_0) { //Метод Хойна
	vector<double> y1, y2, x1, x2;
	double h0 = (B - A) / n;

	//вычисление для h
	x1 = gen_x(h0);
	y1.push_back(Y_A);
	for (int i = 1; i < x1.size(); i++) y1.push_back(y1[i - 1] + (h0 / 2) *
		(function(x1[i - 1], y1[i - 1]) + function(x1[i], y1[i - 1] + h0 * function(x1[i - 1], y1[i - 1]))));
	//----------------

	while (true) {
		n *= 2;
		h0 = (B - A) / n;

		//вычисление для h/2
		x2 = gen_x(h0);
		y2.push_back(Y_A);
		for (int i = 1; i < x2.size(); i++) y2.push_back(y2[i - 1] + (h0 / 2) *
			(function(x2[i - 1], y2[i - 1]) + function(x2[i], y2[i - 1] + h0 * function(x2[i - 1], y2[i - 1]))));
		//------------------

		double acc = accuracy(2, y1, y2); //порядок точности - 2

		if (eps > acc) { //проверка на точность
			*h_0 = h0; //чтобы запомнить точность
			return y2;
		}

		y1.clear();
		x1.clear();

		x1 = x2;
		y1 = y2;

		x2.clear();
		y2.clear();
	}
}

vector<double> Method_advanced_Euler_Cauchy(double n, double eps, double* h_0) { //усовершенствованный метод Эйлера-Коши с итерационной обработкой
	int K = 4;
	double h0 = (B - A) / n;

	vector<double> y1, y2, x1, x2;
	vector<double> y11, y12, y13, y14, y21, y22, y23, y24;

	//вычисление для h
	x1 = gen_x(h0);

	y1.push_back(Y_A);
	for (int i = 0; i < x1.size() - 1; i++) {
		y11.push_back(y1[i] + h0 * function(x1[i], y1[i]));
		y12.push_back(y1[i] + (h0 / 2) * (function(x1[i], y1[i]) + function(x1[i + 1], y11[i])));
		y13.push_back(y1[i] + (h0 / 2) * (function(x1[i], y1[i]) + function(x1[i + 1], y12[i])));
		y14.push_back(y1[i] + (h0 / 2) * (function(x1[i], y1[i]) + function(x1[i + 1], y13[i])));
		y1.push_back(y14[i]);
	}
	//----------------

	while (true) {
		n *= 2;
		h0 = (B - A) / n;

		//вычисление для h/2
		x2 = gen_x(h0);
		y2.push_back(Y_A);
		for (int i = 0; i < x2.size() - 1; i++) {
			y21.push_back(y2[i] + h0 * function(x2[i], y2[i]));
			y22.push_back(y2[i] + (h0 / 2) * (function(x2[i], y2[i]) + function(x2[i + 1], y21[i])));
			y23.push_back(y2[i] + (h0 / 2) * (function(x2[i], y2[i]) + function(x2[i + 1], y22[i])));
			y24.push_back(y2[i] + (h0 / 2) * (function(x2[i], y2[i]) + function(x2[i + 1], y23[i])));
			y2.push_back(y24[i]);
		}
		//------------------

		double acc = accuracy(2, y1, y2); //порядок точности - 2

		if (eps > acc) { //проверка на точность
			*h_0 = h0; //чтобы запомнить точность
			return y2;
		}

		y1.clear();
		x1.clear();
		y11.clear();
		y12.clear();
		y13.clear();
		y14.clear();

		x1 = x2;
		y1 = y2;

		x2.clear();
		y2.clear();
		y21.clear();
		y22.clear();
		y23.clear();
		y24.clear();
	}
}

vector<double> Method_updated_Euler(double n, double eps, double* h_0) { //уточнённый метод Эйлера
	vector<double> y1, y2, x1, x2;
	double h0 = (B - A) / n;

	//вычисление для h
	x1 = gen_x(h0);
	y1.push_back(Y_A); //y[0]
	y1.push_back(y1[0] + h0 * function(x1[0], y1[0])); //по методу эйлера находится y[1]
	for (int i = 1; i < x1.size() - 1; i++) y1.push_back(y1[i - 1] + (2 * h0) * function(x1[i], y1[i]));
	//----------------

	while (true) {
		n *= 2;
		h0 = (B - A) / n;

		//вычисление для h/2
		x2 = gen_x(h0);
		y2.push_back(Y_A);
		y2.push_back(y2[0] + h0 * function(x2[0], y2[0]));
		for (int i = 1; i < x2.size() - 1; i++) y2.push_back(y2[i - 1] + (2 * h0) * function(x2[i], y2[i]));
		//------------------

		double acc = accuracy(2, y1, y2); //порядок точности - 2

		if (eps > acc) { //проверка на точность
			*h_0 = h0; //чтобы запомнить точность
			return y2;
		}
		
		x1.clear();
		y1.clear();

		y1 = y2;
		x1 = x2;

		x2.clear();
		y2.clear();
	}
}

vector<double> Method_mean_point(double n, double eps, double* h_0) { //метод средней точки
	vector<double> y1, y2, x1, x2;
	double h0 = (B - A) / n;

	//вычисление для h
	x1 = gen_x(h0);
	y1.push_back(Y_A); //y[0]
	for (int i = 1; i < x1.size(); i++) y1.push_back(y1[i - 1] + h0 * function(x1[i - 1] + (h0 / 2), y1[i - 1] + (h0 / 2) * function(x1[i - 1], y1[i - 1])));
	//----------------

	while (true) {
		n *= 2;
		h0 = (B - A) / n;

		//вычисление для h/2
		x2 = gen_x(h0);
		y2.push_back(Y_A);
		for (int i = 1; i < x2.size(); i++) y2.push_back(y2[i - 1] + h0 * function(x2[i - 1] + (h0 / 2), y2[i - 1] + (h0 / 2) * function(x2[i - 1], y2[i - 1])));
		//------------------

		double acc = accuracy(2, y1, y2); //порядок точности - 2

		if (eps > acc) { //проверка на точность
			*h_0 = h0; //чтобы запомнить точность
			return y2;
		}

		x1.clear();
		y1.clear();

		y1 = y2;
		x1 = x2;

		x2.clear();
		y2.clear();
	}
}

vector<double> Method_4_order_Runge_Kutta(double n, double eps, double* h_0) { //метод Рунге-Кутты 4-го порядка
	vector<double> eta1, eta2, eta3, eta4;
	vector<double> delta_y;
	double h0 = (B - A) / n;

	vector<double> y1, y2, x1, x2;

	//вычисление для h
	x1 = gen_x(h0);
	y1.push_back(Y_A);

	for (int i = 0; i < x1.size() - 1; i++) {
		eta1.push_back(function(x1[i], y1[i]));
		eta2.push_back(function(x1[i] + (h0 / 2), y1[i] + (h0 / 2) * eta1[i]));
		eta3.push_back(function(x1[i] + (h0 / 2), y1[i] + (h0 / 2) * eta2[i]));
		eta4.push_back(function(x1[i] + h0, y1[i] + h0 * eta3[i]));
		delta_y.push_back((h0 / 6) * (eta1[i] + 2 * eta2[i] + 2 * eta3[i] + eta4[i]));
		y1.push_back(y1[i] + delta_y[i]);
	}
	eta1.clear();
	eta2.clear();
	eta3.clear();
	eta4.clear();
	delta_y.clear();
	//----------------

	while (true) {
		n *= 2;
		h0 = (B - A) / n;

		//вычисление для h/2
		x2 = gen_x(h0);
		y2.push_back(Y_A);

		for (int i = 0; i < x2.size() - 1; i++) {
			eta1.push_back(function(x2[i], y2[i]));
			eta2.push_back(function(x2[i] + (h0 / 2), y2[i] + (h0 / 2) * eta1[i]));
			eta3.push_back(function(x2[i] + (h0 / 2), y2[i] + (h0 / 2) * eta2[i]));
			eta4.push_back(function(x2[i] + h0, y2[i] + h0 * eta3[i]));
			delta_y.push_back((h0 / 6) * (eta1[i] + 2 * eta2[i] + 2 * eta3[i] + eta4[i]));
			y2.push_back(y2[i] + delta_y[i]);
		}
		//------------------

		double acc = accuracy(4, y1, y2); //порядок точности - 4

		if (eps > acc) { //проверка на точность
			*h_0 = h0; //чтобы запомнить точность
			return y2;
		}

		y1.clear();
		x1.clear();

		x1 = x2;
		y1 = y2;

		x2.clear();
		y2.clear();
		eta1.clear();
		eta2.clear();
		eta3.clear();
		eta4.clear();
		delta_y.clear();
	}
}
//--------------------------------------
//вывод в файл таблиц
void output_file(string name_file, vector<double> xx, vector<double> yy, double h_0, double time, string name) {
	//это чтобы для каждого метода создавался свой файл
	//всего будет 6 файлов!!! с таблицами
	ofstream file;
	file.open(name_file);

	file << "Численные методы решения ОДУ первого порядка." << endl;
	file << "Вариант 2.\nЛабораторная работа №3.\ny' = f(x, y) = 2 * (2 - x) * y + 0.01 * exp(-x * x)\nx = [1, 6]\ny(А) = 10" << endl << endl;
	
	file << name << " Время работы этого метода: " << time << endl << endl;
	file << "x\t\ty" << endl;
	for (int i = 0; i < xx.size(); i++) {
		file << xx[i] << "\t\t" << yy[i] << endl;
	}

	file.close();
}

//вывод в файл значений эпсилон и временни для каждого метода
//в файле будет 6 таблиц
void output_(vector<double> epsilon1, vector<double> time1, vector<double> epsilon2, vector<double> time2,
	vector<double> epsilon3, vector<double> time3, vector<double> epsilon4, vector<double> time4,
	vector<double> epsilon5, vector<double> time5, vector<double> epsilon6, vector<double> time6, int n) {
	
	ofstream file;
	file.open("output.txt");

	file << "После проведения " << n << " итераций было замерено время для работы каждого алгоритма \n и учтена точность, введенная пользователем." << endl << endl;

	file << "Для метода Эйлера:" << endl << "epsilon\ttime" << endl;
	for (int i = 0; i < epsilon1.size(); i++) {
		file << epsilon1[i] << "\t" << time1[i] << endl;
	}

	file << "\nДля метода Хойна:" << endl << "epsilon\ttime" << endl;
	for (int i = 0; i < epsilon2.size(); i++) {
		file << epsilon2[i] << "\t" << time2[i] << endl;
	}

	file << "\nДля метода улучшенного Эйлера-Коши:" << endl << "epsilon\ttime" << endl;
	for (int i = 0; i < epsilon3.size(); i++) {
		file << epsilon3[i] << "\t" << time3[i] << endl;
	}

	file << "\nДля метода уточненного Эйлера:" << endl << "epsilon\ttime" << endl;
	for (int i = 0; i < epsilon4.size(); i++) {
		file << epsilon4[i] << "\t" << time4[i] << endl;
	}

	file << "\nДля метода средней точки:" << endl << "epsilon\ttime" << endl;
	for (int i = 0; i < epsilon5.size(); i++) {
		file << epsilon5[i] << "\t" << time5[i] << endl;
	}

	file << "\nДля метода Рунге-Кутты 4-го порядка:" << endl << "epsilon\ttime" << endl;
	for (int i = 0; i < epsilon6.size(); i++) {
		file << epsilon6[i] << "\t" << time6[i] << endl;
	}

	file.close();
}

//вызов всех методов
void calculating(double n, double eps, vector<double> *epsilon1, vector<double>* time1, vector<double>* epsilon2, vector<double>* time2,
				 vector<double>* epsilon3, vector<double>* time3, vector<double>* epsilon4, vector<double>* time4,
				 vector<double>* epsilon5, vector<double>* time5, vector<double>* epsilon6, vector<double>* time6) {
	
	double h0_1, h0_2, h0_3, h0_4, h0_5, h0_6;

	//Эйлер
	double s_time1 = clock();
	vector<double> y1 = Method_Euler(n, eps, &h0_1);
	vector<double> x1 = gen_x(h0_1);
	double e_time1 = clock();
	double second1 = e_time1 - s_time1;
	string f1 = "Метод_Эйлера.txt", name1 = "По методу Эйлера.";
	output_file(f1, x1, y1, h0_1, second1, name1);
	(*epsilon1).push_back(eps);
	(*time1).push_back(second1);

	//Хойн
	double s_time2 = clock();
	vector<double> y2 = Method_Hoina(n, eps, &h0_2);
	vector<double> x2 = gen_x(h0_2);
	double e_time2 = clock();
	double second2 = e_time2- s_time2;
	string f2 = "Метод_Хойна.txt", name2 = "По методу Хойна.";
	output_file(f2, x2, y2, h0_2, second2, name2);
	(*epsilon2).push_back(eps);
	(*time2).push_back(second2);

	//Эйлер-Коши
	double s_time3 = clock();
	vector<double> y3 = Method_advanced_Euler_Cauchy(n, eps, &h0_3);
	vector<double> x3 = gen_x(h0_3);
	double e_time3 = clock();
	double second3 = e_time3 - s_time3;
	string f3 = "Метод_Эйлера_Коши.txt", name3 = "По методу улучшенного Эйлера-Коши.";
	output_file(f3, x3, y3, h0_3, second3, name3);
	(*epsilon3).push_back(eps);
	(*time3).push_back(second3);

	//Уточнённый Эйлер
	double s_time4 = clock();
	vector<double> y4 = Method_updated_Euler(n, eps, &h0_4);
	vector<double> x4 = gen_x(h0_4);
	double e_time4 = clock();
	double second4 = e_time4 - s_time4;
	string f4 = "Метод_Эйлера_уточненный.txt", name4 = "По уточненному методу Эйлера.";
	output_file(f4, x4, y4, h0_4, second4, name4);
	(*epsilon4).push_back(eps);
	(*time4).push_back(second4);

	//Средняя точка
	double s_time5 = clock();
	vector<double> y5 = Method_mean_point(n, eps, &h0_5);
	vector<double> x5 = gen_x(h0_5);
	double e_time5 = clock();
	double second5 = e_time5 - s_time5;
	string f5 = "Метод_средней_точки.txt", name5 = "По методу средней точки.";
	output_file(f5, x5, y5, h0_5, second5, name5);
	(*epsilon5).push_back(eps);
	(*time5).push_back(second5);

	//Рунге-Кутта
	double s_time6 = clock();
	vector<double> y6 = Method_4_order_Runge_Kutta(n, eps, &h0_6);
	vector<double> x6 = gen_x(h0_6);
	double e_time6 = clock();
	double second6 = e_time6 - s_time6;
	string f6 = "Метод_Рунге_Кутты.txt", name6 = "По методу Рунге-Кутты 4-го порядка.";
	output_file(f6, x6, y6, h0_6, second6, name6);
	(*epsilon6).push_back(eps);
	(*time6).push_back(second6);
}


int main() {
	setlocale(LC_ALL, "Russian");

	vector<double> epsilon1, time1, epsilon2, time2, epsilon3, time3, epsilon4, time4, epsilon5, time5, epsilon6, time6;
	int sum = 0;
	char c;

	//создать условие нажатия на клавищу ESC, чтобы сохранять значения эпсилон
	do {
		double eps;
		double n;

		cout << "Введите количество разбиений n: ";
		cin >> n;
		cout << "\nВведите значение epsilon: ";
		cin >> eps;

		calculating(n, eps, &epsilon1, &time1, &epsilon2, &time2, &epsilon3, &time3, &epsilon4, &time4, &epsilon5, &time5, &epsilon6, &time6);

		cout << "Для продолжения работы программы нажмите любую клавишу..." << endl;
		cout << "Для выхода нажмите ESC." << endl;
		c = _getch();

		sum++;
	} while (c != 27);

	output_(epsilon1, time1, epsilon2, time2, epsilon3, time3, epsilon4, time4, epsilon5, time5, epsilon6, time6, sum);

	cout << endl << "Таблицы с рассчитанными координатами x и y сформированы в 6 файлов, которые названы соответствующим образом.\n";
	cout << "Таблицы с различными значениями epsilon выведены в отдельный файл." << endl;

	return 0;
}