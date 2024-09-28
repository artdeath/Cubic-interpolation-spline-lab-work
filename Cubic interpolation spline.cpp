#include <iostream>
#include <vector>
#include <algorithm>


class Spline {

private:

	double r;

	int n; //количество элементов (количество сегментов - 1)

	std::vector<double> h;

	std::vector<std::pair<double, double>> Grid;

	std::vector<double> a; std::vector<double> b; std::vector<double> c; std::vector<double> d; //коэффициенты сплайна

public:

	//конструктор; для регулярной сетки r = 1. Отрезок задаётся начальной точкой, шагом и параметром релаксации
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

		return std::pow(x, 3) + x;

	}

	//значение производной аналитической функции
	double df(double x) {

		return 3 * std::pow(x, 2) + 1.;

	}

	//значение второй производной аналитической функции

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

	//значение его производной
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

	//значение его второй производной

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

	//вычисляет коэффициенты сплайна
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

		double y; //знаменатель метода прогонки, для удобства

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

			std::reverse(c.begin(), c.end()); //~~вычисляет значения с конца

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

			c.push_back(0);
			
			std::reverse(c.begin(), c.end());

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

	Spline X(0.125, -1., 1., 9);

	X.gMake();

	double r;

	int k = 0;

	std::vector<double> P;

	P.push_back(-1.); P.push_back(0.); P.push_back(-0.5);

	while (k < 8) {

		bool j = true;

		do {

			r = static_cast<double>(rand()) / RAND_MAX;

		} while (r == 0. || r == 1.);

		r = r * std::fabs(X.get_xGrid_at(0) - X.get_xGrid_at(X.get_n() - 1)) + X.get_xGrid_at(0); //сдвиг диапазона

		for (int i = 0; i < X.get_n(); i++) {

			if (std::fabs(r - X.get_xGrid_at(i)) < std::numeric_limits<double>::epsilon()) {

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

	//Вывод сетки

	std::cout << std::scientific << "h = 0.125" << std::endl;

	for (int i = 0; i < X.get_n(); i++) {


		std::cout << std::fixed << X.get_xGrid_at(i) << '\t' << std::scientific << X.get_fGrid_at(i) << std::endl;

	}

	std::cout << std::endl;

	//Вывод функции

	std::cout << std::scientific << "x" << "\t\t" << "f(x)" << "\t\t" << "f'(x)" << "\t\t" << "f''(x)" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::fixed << P.at(i) << std::scientific << '\t' << X.f(P.at(i)) << '\t' << X.df(P.at(i)) << '\t' << X.ddf(P.at(i)) << std::endl;

	}

	std::cout << std::endl;

	//Вывод сплайна

	std::cout << std::scientific << "x" << "\t\t" << "g(x)" << "\t\t" << "g'(x)" << "\t\t" << "g''(x)" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::fixed << P.at(i) << std::scientific << '\t' << X.g(P.at(i)) << '\t' << X.dg(P.at(i)) << '\t' << X.ddg(P.at(i)) << std::endl;

	}

	std::cout << std::endl;

	//Погрешности

	std::cout << std::scientific << "h = 0.125" << std::endl;
	std::cout << std::scientific << "|f - g|\t\t|f' - g'|" << std::endl;

	for (int i = 0; i < 11; i++) {

		std::cout << std::scientific << std::fabs(X.g(P.at(i)) - X.f(P.at(i))) << '\t' << std::fabs(X.dg(P.at(i)) - X.df(P.at(i))) << std::endl;

	}

	std::cout << std::endl;


}
