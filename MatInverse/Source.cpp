#include "Twin.h"
#include <sstream>

std::vector<std::vector<double>> getA(Twin&);
void dmat(Twin&, const std::vector<std::vector<double>>&);
std::vector<std::vector<double>> matinv(const std::vector<std::vector<double>>&);

int main() {

	Twin t("Matrix Inverse");

	while (true) {

		std::vector<std::vector<double>> A = getA(t);

		std::vector<std::vector<double>> Ainv = matinv(A);
		t.println("Inverse of Matrix A:");
		t.println();

		dmat(t, Ainv);
		t.println("Continue? [y/n]");
		t.display();

		std::string input;
		std::getline(std::cin, input);

		if (input != "y" && input != "Y") break;
	}

	return EXIT_SUCCESS;
}

std::vector<std::vector<double>> getA(Twin& t) {

	std::vector<std::vector<double>> A;

	bool enter = true;
	while (enter) {

		while (true) {
			t.println("Enter the first row of matrix A");
			t.println("A must be a square, invertible matrix.");
			t.println("Example: 1 2 3 4");
			t.display();

			double in;
			std::string input;
			std::getline(std::cin, input);

			if (input == "") {
				t.println("Try again");
				t.println();
			}
			else {
				std::stringstream iss(input);
				A.push_back(std::vector<double>());
				while (iss >> in) A[0].push_back(in);
				break;
			}
		}

		int m = 1;

		while (m < int(A[0].size())) {

			dmat(t, A);

			t.println("Enter the next row of the matrix.");
			t.display();

			double in;
			std::string input;
			std::getline(std::cin, input);
			std::stringstream iss(input);

			A.push_back(std::vector<double>());
			while (iss >> in) A[m].push_back(in);

			if (A[m].size() == A[0].size())
				m++;
			else {
				A.pop_back();
				t.println("Invalid row.");
				t.println();
			}
		}

		t.println("Matrix A:");
		t.println();
		dmat(t, A);
		t.println("Save and continue? [y/n]");
		t.display();

		std::string input;
		std::getline(std::cin, input);

		if (input != "y" && input != "Y") {
			A.clear();
		}
		else enter = false;
	}
	return A;
}

void dmat(Twin& t, const std::vector<std::vector<double>>& mat) {

	t.println();

	for (int i = 0; i < int(mat.size()); i++) {

		t.print("[ ");

		for (int j = 0; j < int(mat[0].size()); j++) {
			t.print(mat[i][j]);
			if (j < int(mat[0].size()) - 1) t.print(",");
			t.print(" ");
		}

		t.print("]");
		t.println();

	}
}

std::vector<std::vector<double>> matinv(const std::vector<std::vector<double>>& A) {

	int n = A.size();

	std::vector<std::vector<double>> l, u(A), Ainv(n);

	for (int i = 0; i < n; i++) {
		l.push_back(std::vector<double>(n, 0));
		l[i][i] = 1;
	}

	for (int i = 0; i < n; i++) {

		for (int k = i + 1; k < n; k++) {
			l[k][i] = u[k][i] / u[i][i];

			for (int j = i; j < n; j++) u[k][j] += -l[k][i] * u[i][j];

		}
	}

	for (int j = 0; j < n; j++) {

		std::vector<double> x(n), y(n), b(n,0);

		b[j] = 1;

		for (int i = 0; i < n; i++) {
			y[i] = b[i] / l[i][i];
			for (int k = 0; k < n; k++) {
				b[k] -= l[k][i] * y[i];
			}
		}

		for (int i = n - 1; i >= 0; i--) {
			x[i] = y[i] / u[i][i];
			for (int k = i - 1; k >= 0; k--) {
				y[k] -= u[k][i] * x[i];
			}
		}

		for (int k = 0; k < n; k++) 
			Ainv[k].push_back(x[k]);
	}

	return Ainv;

}