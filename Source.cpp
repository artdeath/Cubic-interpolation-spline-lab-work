#include <iostream>
#include <vector>
#include <algorithm>


class Spline {

private:

	double r;

	int n; //количество сегментов (количество элементов - 1)

	std::vector<double> h;

	std::vector<std::pair<double, double>> Grid;

	std::vector<double> a; std::vector<double> b; std::vector<double> c; std::vector<double> d; //массивы коэффициентов сплайна

public:

	//конструктор для двух сегментов из задания, default)
	Spline() {

		h.push_back(1.);

		n = 2;

		double start = -1.;

		for (int i = 0; i < 3; i++) {

			Grid.push_back(std::pair<double, double>(start, f(start)));

			start += h.back();

		}


	};

	//конструктор с параметрами; для регулярной сетки r = 1
	Spline(double input_h, double input_start, double input_r, int input_n) {

		h.push_back(input_h);

		r = input_r;

		n = input_n;

		double start = input_start;

		if (std::abs(1. - input_r) < std::numeric_limits<double>::epsilon()) {

			for (int i = 0; i < input_n; i++) {

				Grid.push_back(std::pair<double, double>(start, f(start)));

				start += h.back();

			}

		}
		else {

			Grid.push_back(std::pair<double, double>(start, f(start)));

			for (int i = 0; i < input_n - 1; i++) {

				start += h.back();

				Grid.push_back(std::pair<double, double>(start, f(start)));

				h.push_back(h.back() * r);

			}

			h.pop_back();

		}


	};

	//возвращает значение аналитической функции
	double f(double x) {

		return std::pow(x, 3) - x;

		//return std::exp(x);

	}

	//возвращает значение первой производной
	double df(double x) {

		return 3 * std::pow(x, 2) - 1.;

	}

	//возвращает значение второй производной

	double ddf(double x) {

		return 6 * x;

	}

	//возвращает значение сплайна
	double g(double x) {

		double Grid_Inc = Grid.at(0).first + h.at(0);

		if (std::abs(1. - r) >= std::numeric_limits<double>::epsilon()) {

			for (int i = 0; i < n - 1; i++) {

				if ((x >= Grid.at(0).first) and (x <= Grid_Inc)) {

					return a.at(i) + b.at(i) * (x - Grid.at(i).first) + c.at(i) * std::pow((x - Grid.at(i).first), 2) + d.at(i) * std::pow((x - Grid.at(i).first), 3);

				}

				Grid_Inc += h.at(i);

			}
		}
		else {

			for (int i = 0; i < n - 1; i++) {
				
				if ((x >= Grid.at(0).first) and (x <= Grid_Inc)) {

					return a.at(i) + b.at(i) * (x - Grid.at(i).first) + c.at(i) * std::pow((x - Grid.at(i).first), 2) + d.at(i) * std::pow((x - Grid.at(i).first), 3);

				}

				Grid_Inc += h.back();

			}

		}

		return 0;

	}

	//возвращает значение первой производной
	double dg(double x) {

		double Grid_Inc = Grid.at(0).first + h.at(0);

		if (std::abs(1. - r) >= std::numeric_limits<double>::epsilon()) {

			for (int i = 0; i < n - 1; i++) {

				if ((x >= Grid.at(0).first) and (x <= Grid_Inc)) {

					return 3 * d.at(i) * std::pow(x, 2) + (2 * c.at(i) - 6 * d.at(i) * Grid.at(i).first) * x + 3 * d.at(i) * std::pow(Grid.at(i).first, 2) - 2 * c.at(i) * Grid.at(i).first + b.at(i);

				}

				Grid_Inc += h.at(i);

			}
		}
		else {

			for (int i = 0; i < n - 1; i++) {

				if ((x >= Grid.at(0).first) and (x <= Grid_Inc)) {

					return 3 * d.at(i) * std::pow(x, 2) + (2 * c.at(i) - 6 * d.at(i) * Grid.at(i).first) * x + 3 * d.at(i) * std::pow(Grid.at(i).first, 2) - 2 * c.at(i) * Grid.at(i).first + b.at(i);

				}

				Grid_Inc += h.back();

			}

		}

		return 0;

	}

	//возвращает значение второй производной

	double ddg(double x) {

		double Grid_Inc = Grid.at(0).first + h.at(0);

		if (std::abs(1. - r) >= std::numeric_limits<double>::epsilon()) {

			for (int i = 0; i < n - 1; i++) {

				if ((x >= Grid.at(0).first) and (x <= Grid_Inc)) {

					return 6 * d.at(i) * (x - Grid.at(i).first) + 2 * c.at(i);

				}

				Grid_Inc += h.at(i);

			}
		}
		else {

			for (int i = 0; i < n - 1; i++) {

				if ((x >= Grid.at(0).first) and (x <= Grid_Inc)) {

					return 6 * d.at(i) * (x - Grid.at(i).first) + 2 * c.at(i);

				}

				Grid_Inc += h.back();

			}

		}

		return 0;


	}

	//вычисляет векторы a, b, c, d
	void gMake() {
		
		std::vector<double> B;

		std::vector<double> F;

		if (std::abs(1. - r) < std::numeric_limits<double>::epsilon()) {

			for (int i = 0; i < n - 2; i++) {

				B.push_back(4 * h.back());

				F.push_back(3 * (((Grid.at(i + 2).second - Grid.at(i + 1).second) / h.back()) - ((Grid.at(i + 1).second - Grid.at(i).second) / h.back())));


			}

		}
		else {

			for (int i = 0; i < n - 2; i++) {

				B.push_back(2 * (h.at(i) + h.at(i + 1)));

				F.push_back(3 * (((Grid.at(i + 2).second - Grid.at(i + 1).second) / h.at(i + 1)) - ((Grid.at(i + 1).second - Grid.at(i).second) / h.at(i))));

			}


		}

		std::vector<double> v; std::vector<double> u;

		double y; //знаменатель

		if (std::abs(1. - r) >= std::numeric_limits<double>::epsilon()) {

			y = B.at(0); v.push_back(-h.at(1) / y); u.push_back(F.at(0) / y);

			if (n > 3) {

				for (int i = 1; i < n - 3; i++) {

					y = B.at(i) + h.at(i) * v.at(i - 1);

					v.push_back(-h.at(i + 1) / y);

					u.push_back((F.at(i) - h.at(i) * u.at(i - 1)) / y);

				}

				y = B.at(n - 3) + h.at(n - 3) * v.at(n - 4); v.push_back(0.); u.push_back(F.at(n - 3) - h.at(n - 3) * u.at(n - 4));

			}

			c.push_back(u.at(n - 3));

			if (n > 3) {

				for (int i = n - 4; i >= 0; i--) {

					c.push_back(v.at(i) * c.back() + u.at(i));

				}

			}

			c.push_back(0.); //метод прогонки~~

			std::reverse(c.begin(), c.end()); //~~записывает значение в обратном порядке

			for (int i = 0; i < n - 2; i++) {

				a.push_back(Grid.at(i).second);

				b.push_back(((Grid.at(i + 1).second - Grid.at(i).second) / h.at(i)) - ((2 * c.at(i) + c.at(i + 1)) * h.at(i) / 3));

				d.push_back((c.at(i + 1) - c.at(i)) / (3 * h.at(i)));

			}

			a.push_back(Grid.at(n - 2).second);

			b.push_back(((Grid.at(n - 1).second - Grid.at(n - 2).second) / h.at(n - 2)) - ((2 * c.at(n - 2)) * h.at(n - 2) / 3));

			d.push_back(-c.at(n - 2) / (3 * h.at(n - 2)));

		}
		else {

			y = B.at(0); v.push_back(-h.back() / y); u.push_back(F.at(0) / y);

			if (n > 3) {

				for (int i = 1; i < n - 3; i++) {

					y = B.at(i) + h.back() * v.at(i - 1);

					v.push_back(-h.back() / y);

					u.push_back((F.at(i) - h.back() * u.at(i - 1)) / y);

				}

				y = B.at(n - 3) + h.back() * v.at(n - 4); v.push_back(0.); u.push_back(F.at(n - 3) - h.back() * u.at(n - 4));

			}

			c.push_back(u.at(n - 3));

			if (n > 3) {

				for (int i = n - 4; i >= 0; i--) {

					c.push_back(v.at(i) * c.back() + u.at(i));

				}

			}

			c.push_back(0); //метод прогонки~~

			std::reverse(c.begin(), c.end()); //~~записывает значение в обратном порядке

			for (int i = 0; i < n - 2; i++) {

				a.push_back(Grid.at(i).second);

				b.push_back(((Grid.at(i + 1).second - Grid.at(i).second) / h.back()) - ((2 * c.at(i) + c.at(i + 1)) * h.back() / 3));

				d.push_back((c.at(i + 1) - c.at(i)) / (3 * h.back()));

			}

			a.push_back(Grid.at(n - 2).second);

			b.push_back(((Grid.at(n - 1).second - Grid.at(n - 2).second) / h.back()) - ((2 * c.at(n - 2)) * h.back() / 3));

			d.push_back(-c.at(n - 2) / (3 * h.back()));

		}

	}

	double get_n() {

		return n;

	}

	double get_xGrid_at(int i) {

		return Grid.at(i).first;

	}

	double get_fGrid_at(int i) {

		return Grid.at(i).second;

	}

};

int main() {

	std::srand(static_cast<unsigned int>(time(0)));

	Spline X1(0.5, -1., 1., 3);

	X1.gMake();

	Spline X2(0.25, -1., 1., 5);

	X2.gMake();

	Spline X3(0.125, -1., 1., 9);

	X3.gMake();

	double r;

	int k = 0;

	std::vector<double> P;

	P.push_back(-1.); P.push_back(0.); P.push_back(-0.5);

	while (k < 8) {

		bool j = true;

		do {

			r = static_cast<double>(rand()) / RAND_MAX;

		} while (r == 0. || r == 1.);

		r -= 1; //сдвигаме диапазон

		//r /= 10.;

		for (int i = 0; i < X3.get_n(); i++) {

			if (r == X3.get_xGrid_at(i)) {

				j = false;

				break;

			}

		}

		if (j) {

			k++;

			P.push_back(r);

		}

	}

	std::sort(P.begin(), P.end());

	//вывод сеток

	std::cout << std::scientific << "h = 0.5\t\t\t\th = 0.25\t\t\th = 0.125" << std::endl;
	
	for (int i = 0; i < X3.get_n(); i++) {

		if (i < X1.get_n()) {

			std::cout << std::fixed << X1.get_xGrid_at(i) << '\t' << std::scientific << X1.get_fGrid_at(i) << '\t';

		}
		else {

			std::cout << "\t\t\t\t";

		}

		if (i < X2.get_n()) {

			std::cout << std::fixed << X2.get_xGrid_at(i) << '\t' << std::scientific << X2.get_fGrid_at(i) << '\t';

		}
		else {

			std::cout << "\t\t\t\t";

		}

		std::cout << std::fixed << X3.get_xGrid_at(i) << '\t' << std::scientific << X3.get_fGrid_at(i) << std::endl;

	}

	std::cout << std::endl;

	//вывод функции

	std::cout << std::scientific << "x" << "\t\t" << "f(x)" << "\t\t" << "f'(x)" << "\t\t" << "f''(x)" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::fixed << P.at(i) << std::scientific << '\t' << X1.f(P.at(i)) << '\t' << X1.df(P.at(i)) << '\t' << X1.ddf(P.at(i)) << std::endl;

	}

	std::cout << std::endl;

	//первое разбиение

	std::cout << std::scientific << "x" << "\t\t" << "g(x)" << "\t\t" << "g'(x)" << "\t\t" << "g''(x)" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::fixed << P.at(i) << std::scientific << '\t' << X1.g(P.at(i)) << '\t' << X1.dg(P.at(i)) << '\t' << X1.ddg(P.at(i)) << std::endl;

	}

	std::cout << std::endl;

	//второе разбиение

	std::cout << std::scientific << "x" << "\t\t" << "g(x)" << "\t\t" << "g'(x)" << "\t\t" << "g''(x)" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::fixed << P.at(i) << std::scientific << '\t' << X2.g(P.at(i)) << '\t' << X2.dg(P.at(i)) << '\t' << X2.ddg(P.at(i)) << std::endl;

	}

	std::cout << std::endl;

	//третье разбение

	std::cout << std::scientific << "x" << "\t\t" << "g(x)" << "\t\t" << "g'(x)" << "\t\t" << "g''(x)" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::fixed << P.at(i) << std::scientific << '\t' << X3.g(P.at(i)) << '\t' << X3.dg(P.at(i)) << '\t' << X3.ddg(P.at(i)) << std::endl;

	}

	std::cout << std::endl;

	//Погрешности

	std::cout << std::scientific << "h = 0.5\t\t\t\th = 0.25\t\t\th = 0.125" << std::endl;
	std::cout << std::scientific << "|f - g|\t\t|f' - g'|\t|f - g|\t\t|f' - g'|\t|f - g|\t\t|f' - g'|" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::scientific << std::fabs(X1.g(P.at(i)) - X1.f(P.at(i))) << '\t' << std::fabs(X1.dg(P.at(i)) - X1.df(P.at(i))) << '\t';
		std::cout << std::scientific << std::fabs(X2.g(P.at(i)) - X2.f(P.at(i))) << '\t' << std::fabs(X2.dg(P.at(i)) - X2.df(P.at(i))) << '\t';
		std::cout << std::scientific << std::fabs(X3.g(P.at(i)) - X3.f(P.at(i))) << '\t' << std::fabs(X3.dg(P.at(i)) - X3.df(P.at(i))) << std::endl;

	}

	std::cout << std::endl;


}