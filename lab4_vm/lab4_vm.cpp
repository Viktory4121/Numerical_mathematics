// Вариант 2.

#include <iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cstdlib>
#include<string>

using namespace std;

//---------------------ПО ВАРИАНТУ---------------------
const double A = 1.0;
//инициализация границ для x и t происходит в функции calculating
//---------------------------
//f(t)
double f_t(double x, double t, double coef) { 
	return coef * (5 / (2 * pow(x + 1.5 * t + 1, 3)));
}

//u(0,x)
vector<double> u_0_x(vector<double> x) {
	vector<double> u;
	for (int i = 0; i < x.size(); i++) {
		u.push_back(1 / (1 + x[i]));
	}
	return u;
}

//du/dt(0,x)
double du_dt_0_x(double x) { 
	return (-3 / (2 * (1 + x) * (1 + x)));
}

//u(t,0)
vector<double> u_t_0(vector<double> t) { 
	vector<double> u;
	for (int i = 1; i < t.size(); i++) {
		u.push_back(1 / (1.5 * t[i] + 1));
	}
	return u;
}

//u(t,1)
vector<double> u_t_1(vector<double> t) { 
	vector<double> u;
	for (int i = 1; i < t.size(); i++) {
		u.push_back(1 / (1.5 * t[i] + 2));
	}
	return u;
}
//-----------------------------------------------------
//Генерация переменных в соответствие с их границами и заданным шагом
vector<double> gen_var(double step, double a, double b) {
	vector<double> var;

	for (double i = a; i <= b; i += step) {
		var.push_back(i);
	}
	//if (var[var.size() - 1] != b) var.push_back(b);
	return var;
}

//Вывод таблицы в файл
void output_table(vector<vector<double>> u, vector<double> x, vector<double> t, const char* name_file, const char* message) {
	FILE* file;
	fopen_s(&file, name_file, "w");

	fprintf_s(file, "Лабораторная работа №4. Вариант 2.");
	fprintf_s(file, "\nПосле проведения вычисления, получим таблицу:\n");
	fprintf_s(file, message);


	int tt = t.size() - 1;
	for (int i = 0; i < t.size(); i++) {
		fprintf_s(file, "%7.4f|", t[tt]);							//вывод столбца t

		for (int j = 0; j < x.size(); j++) {
			fprintf_s(file, "%11.3f", u[i][j]);						//вывод значений функции u	
		}
		fprintf_s(file, "\n");
		tt--;
	}

	fprintf_s(file, " t / x |");
	for (int i = 0; i < x.size(); i++) {
		fprintf_s(file, "%11.3f", x[i]);								//вывод строки x		
	}
	fclose(file);
}

//вызов основных функций
void calculating(double h, double tau, double T, double coef, bool logic) {
	const double Ax = 0.0;		//границы для х
	const double Bx = 1.0;
	const double At = 0.0;		//границы для t
	const double Bt = T;

	vector<double> x = gen_var(h, Ax, Bx);
	vector<double> t = gen_var(tau, At, Bt);

	//Далее заполняется в 5 этапов таблица
	//сначала 1, 3, 2 (см. на рисунок в блокноте)
	vector<double> table_1 = u_0_x(x);
	vector<double> table_2 = u_t_1(t);
	vector<double> table_3 = u_t_0(t);

	//затем (4) часть в вектор
	vector<double> table_4;
	for (int i = 1; i < x.size() - 1; i++) {
		table_4.push_back(table_1[i] + tau * du_dt_0_x(x[i]) + ((tau * tau) / 2) * (A * A * (2 / (pow(1 + x[i], 3)))) + f_t(x[i], t[0], coef));
	}

	//Создадим окончательную таблицу, в которой будут помещеы ранее вычисленные значения:
	//t сверху вниз идёт от tm до t0 (из конца в начало)
	//х слева направо идёт из начала в конец
	vector<vector<double>> table(t.size(), vector<double>(x.size()));

	int k = t.size() - 2;
	for (int i = 0; i < t.size() - 1; i++) {
		table[i][0] = table_3[k];				//заполнение первого столбца таблицы до t1
		table[i][x.size() - 1] = table_2[k];	//заполнение последнего столбца таблицы до t1
		k--;
	}
	table_2.clear();
	table_3.clear();

	for (int i = 0; i < x.size(); i++) {
		table[t.size() - 1][i] = table_1[i];			//заполнение t0 строки таблицы
	}
	table_1.clear();

	k = 0;
	for (int i = 1; i < x.size() - 1; i++) {
		table[t.size() - 2][i] = table_4[k];	//заполнение t1 строки с х1 элемента до х_n-1
		k++;
	}
	table_4.clear();


	//теперь заполним центр таблицы (5) в общую таблицу
	for (int i = t.size() - 3; i > 0; i--) {
		for (int j = 1; j <= x.size() - 2; j++) {
			table[i][j] = A * A * tau * tau * ((table[i + 1][j + 1] - 2 * table[i + 1][j] + table[i + 1][j - 1]) / (h * h)) 
						+ f_t(x[j], t[i + 1], coef) * tau * tau + 2 * table[i + 1][j] - table[i + 2][j];
		}
	}

	if (coef == 1 && logic) {
		output_table(table, x, t, "case1.txt", "при условии, что выполнено условие Куранта и коэффициент при функции равен 1.\n\n");
	} else if (coef == 1 && !logic){
		output_table(table, x, t, "case2.txt", "при условии, что не выполнено условие Куранта и коэффициент при функции равен 1.\n\n");
	} else if (coef != 1 && logic) {
		output_table(table, x, t, "case3.txt", "при условии, что выполнено условие Куранта и коэффициент при функции не равен 1.\n\n");
	} else {
		output_table(table, x, t, "case4.txt", "при условии, что не выполнено условие Куранта и коэффициент при функции не равен 1.\n\n");
	}
	

	x.clear();
	t.clear();
	table.clear();
}

int main() {
    setlocale(LC_ALL, "Russian");

	double tau;
	double h;
	double T;
	double coef;
	bool logic = false;

	cout << "Введите правую границу T для t: ";
	cin >> T;
	cout << "Введите шаг h для x: ";
	cin >> h;
	cout << "Введите шаг tau для t: ";
	cin >> tau;

	if (A * (tau / h) <= 1) {
		cout << "\nВыполнено условие Куранта.\n";
		logic = true;
	} else {
		cout << "\nУсловие Куранта не выполнено.\n";
	}

	cout << "Введите коэффициент при функции f(t): ";
	cin >> coef;

	calculating(h, tau, T, coef, logic);

	cout << "\nПроизведён вывод в файл с результатом выполнения программы.";

    return 0;
}