#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <conio.h>
#include <string>

double polynomial(double x)
{
    return 0;
}

struct BisectionMethod{  
    double Method(double left, double right, double epsilion)
    {
        double leftValue = polynomial(left);
        double rightValue = polynomial(right);
        double mid = left+abs(left-right) / 2;
        double midValue = polynomial(mid);
        while ( abs(midValue) > epsilion) {
            if ((leftValue * midValue) < 0)
                right = mid;
            else
                left = mid;

            mid = left + abs(left - right) / 2;
            leftValue = polynomial(left);
            rightValue = polynomial(right);
            midValue = polynomial(mid);
        }

        return mid;
    }
    };




struct LU
{
    void croutsMethod(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& U) {
        int n = A.size();

        for (int i = 0; i < n; ++i) {
            U[i][i] = 1;

            for (int j = i; j < n; ++j) {
                double sum = 0;

                for (int k = 0; k < i; ++k) {
                    sum += L[i][k] * U[k][j];
                }

                U[i][j] = A[i][j] - sum;
            }

            for (int j = i; j < n; ++j) {
                if (i == j) {
                    L[i][j] = 1;
                }
                else {
                    double sum = 0;

                    for (int k = 0; k < i; ++k) {
                        sum += L[j][k] * U[k][i];
                    }

                    L[j][i] = (A[j][i] - sum) / U[i][i];
                }
            }
        }
    }

    void solveSystem(std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& U, std::vector<double>& B, std::vector<double>& X) {
        int n = L.size();

        std::vector<double> Y(n, 0);

        // Solve LY = B
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < i; ++j) {
                sum += L[i][j] * Y[j];
            }
            Y[i] = (B[i] - sum) / L[i][i];
        }

        // Solve UX = Y 
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < n; ++j)
            {
                sum += U[i][j] * X[j];
             
             }
            X[i] = (Y[i] - sum) / U[i][i];
        }
    }
};

struct SimpleIterationForMatrix
{
private:
    double expression(std::vector<double> vec, std::vector<double> ficent, double b)
    {
        return vec[0] * ficent[0] + ficent[1] * vec[1] + ficent[2] * vec[2] + b;
    }
public:
    std::vector<double> Method(std::vector<std::vector<double>> matrix, std::vector<double> b, double epsilion)
    {
        std::vector<double> x0 = { 0,0,0 };
        std::vector<double> x1(3);
        for (int i = 0; i < x1.size(); ++i)
        {
            x1[i] = expression(x0, matrix[i], b[i]);
        }
        for (auto& x : x1)
            std::cout << x << " ";
        std::cout << "\n";

        double difx0 = abs(x0[0] - x1[0]);
        double difx1 = abs(x0[1] - x1[1]);
        double difx2 = abs(x0[2] - x1[2]);

        double diff = std::max(std::max(difx0, difx1), difx2);
        while (diff > epsilion)
        {
            for (auto& x : x0)
                std::cout << x << " ";
            std::cout << "\n";
            x0 = x1;
            for (int i = 0; i < x1.size(); ++i)
            {
                x1[i] = expression(x0, matrix[i], b[i]);
            }

            difx0 = abs(x0[0] - x1[0]);
            difx1 = abs(x0[1] - x1[1]);
            difx2 = abs(x0[2] - x1[2]);

            diff = std::max(std::max(difx0, difx1), difx2);
            std::cout << "\ndiff " << diff << "\n";
        }
        return x1;

    }

};


int main()
{



    std::vector<std::vector<double>> a ={{0.1, -0.3,0.1},{0.4, 0.2, -0.3},{0.3,-0.1,0.5}};
    std::vector<double> b = { -1, 2,1};
   
    SimpleIterationForMatrix ex1;
    std::vector<double> result = ex1.Method(a, b, 0.0001);
    for (auto& x : result)
          std::cout << x << " ";


          
    LU lu;
    int n = a.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0));
    std::vector<double> X(n, 0);
    lu.croutsMethod(a, L, U);
    lu.solveSystem(L, U, b, X);
    std::cout << "\nSolution With LU: ";
    for (const auto& x : X)
    {
        if (x == -0) std::cout << abs(x) << " ";
         else std::cout << x << " ";
    }

}
