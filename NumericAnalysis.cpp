#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <conio.h>
#include <string>
double polynomial(double x)
{
    //return (x - 6) * (x - 5) * (x - 4) * (x - 3) * (x - 2) * (x - 1);
    //return x * x * x - x + 1;
    //return x * x * x * x - 2 * x * x * x + 4 * x - 1;
    // return pow(x,4) + 3*pow(x,3) - 2*x;
    //return 6 * pow(x, 5) - 16*pow(x, 3)  - 3;
    /*Simple Iteration*/ //return x - ((x * x * x + 3 * x * x - 7) / 13);
    /*Newton*/// return (x * x * x + 3 * x * x - 7) / (3 * x * x + 6 * x);
    //double y = x - 0.01;
    //return (x * x * x + 3 * x * x - 7) / ((x * x * x + 3 * x * x - 7) - (x * x * x + 3 * y * y - 7));

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

struct SimpleIteration{
//g(x) = x - (f(x)/Alpha), |Alpha| > 1/2(max[a,b]|f'(x)|) 
    double Iteration(double xk_1, double epsilion)
    {
        double xk = polynomial(xk_1);
        while (abs(xk - xk_1) > epsilion)
        {
            xk_1 = xk;
            xk = polynomial(xk);
        }
        return xk;
    }
};

struct NewtonMethod{
private:
    double f2Dif(double x)
    {
        return 6 * x + 6;
    }
    double FFunc(double x)
    {
        return x * x * x + 3 * x * x - 7;
    }
public:
    double Method(double a, double b,double epsilion)
    {
        double xk_1 = a;
        double xk = 99;
        if (f2Dif(a) * FFunc(a)> 0)
        {
            xk_1 = a;
        }
        else
            xk_1 = b;
        while (abs(FFunc(xk)) > epsilion)
        {
            xk = xk_1 - polynomial(xk_1);
            xk_1 = xk;
        }
        return xk;

          
    }
};

struct IntersectionMethod{
private:
    double f2Dif(double x)
    {
        return 6 * x + 6;
    }
    double FFunc(double x)
    {
        return x * x * x + 3 * x * x - 7;
    }
public:
    double method(double a, double b, double epsilion)
    {
        double xk = 0;
        double xk_1 = 199;
        double xk_2 = 299;
        if (f2Dif(a) * FFunc(a) > 0)
        {
            xk_1 = a;
        }
        else
            xk_1 = b;
        if (xk_1 == a)
        {
            xk_2 = a + 0.001;
        }
        else xk_2 = b - 0.001;

            while (abs(FFunc(xk)) > epsilion)
            {
                xk = xk_1 - ((xk_1 * xk_1 * xk_1 + 3 *xk_1 * xk_1 - 7) / ((xk_1 * xk_1 * xk_1 + 3 * xk_1 * xk_1 - 7) - (xk_2 * xk_2 * xk_2 + 3 * xk_2 * xk_2 - 7))*(xk_1-xk_2));
                xk_2 = xk_1;
                xk_1 = xk;
            }
        return xk;
    }

};

struct Gauss {
    std::vector<double> GaussMethod(std::vector<std::vector<double>> A, std::vector<double> b)
    {
        std::vector<double> X(A.size(), 1);

        for (int i = 0; i < A.size(); ++i)
        {
            for (int j = i + 1; j < A.size(); ++j)
            {
                if (A[i][i] == 0) throw "Sorry I can't solve your problem";
                double factor = A[j][i] / A[i][i];
                for (int k = 0; k < A[0].size(); ++k)
                {
                    A[j][k] -= factor * A[i][k];
                }
                b[j] -= factor * b[i];

            }
        }


        for (int i = A.size() - 1; i >= 0; --i)
        {
            double sum = 0;
            for (int j = i + 1; j < A[0].size(); ++j)
            {
                sum += A[i][j] * X[j];
            }
            if (A[i][i] == 0) throw "Sorry I can't solve your problem";
            X[i] = (b[i] - sum) / A[i][i];
        }
        return X;

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
            for (int j = i + 1; j < n; ++j) { sum += U[i][j] * X[j]; }
            X[i] = (Y[i] - sum) / U[i][i];
        }
    }
};

struct SimpleIterationForMatrix {
   
private:
    double expression(std::vector<double> vec, std::vector<double> ficent, double b)
    {
        return vec[0] * ficent[0] + ficent[1] * vec[1] + ficent[2] * vec[2]+b;
    }
public:
    std::vector<double> Method(std::vector<std::vector<double>> matrix,std::vector<double> b, double epsilion )
    {
        std::vector<double> x0 = { 0,0,0 };
        std::vector<double> x1 = { 1,1,1 };
        for(int i = 0; i < x1.size();++i)
        {
            x1[i] = expression(x0, matrix[i], b[i]);
        }

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

        }
        return x1;

    }

};

int main()
{
	std::vector<std::vector<double>> a ={{0.1, -0.3,0.1},{0.4, 0.2, -0.3},{0.3,-0.1,0.5}};
	std::vector<double> b = { -1, 2,1};

	SimpleIterationForMatrix ex1;
	std::vector<double> result = ex1.Method(a, b, 0.01);
    
	for (auto& x : result)
        	std::cout << x << " ";

	
	std::vector<std::vector<double>> a ={{-2, 3, 6},{4, -1, 3},{1, 7, -1}};
	std::vector<double> b = { 5, -5,6};

    // Gauss
    Gauss G;
	std::vector<double> Y = G.GaussMethod(a, b);
    std::cout << "\nSolution with Gauss: ";
	for (const auto& x : Y)
	{
		if (x == -0) std::cout << abs(x) <<" ";
		else std::cout << x <<" ";
	}
    
    // LU
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
    
    //Bisection
    //BisectionMethod BisecMethod;
    //std::cout << "\nSolution Bisection Method For Polynominal: " << BisecMethod.Method(-10  , -2.14, 0.01);
    //std::cout << "\nSolution Bisection Method For Polynominal: " << BisecMethod.Method(-2.14, -0.54, 0.01);
    //std::cout << "\nSolution Bisection Method For Polynominal: " << BisecMethod.Method(-0.54,  0.43, 0.01);
    //std::cout << "\nSolution Bisection Method For Polynominal: " << BisecMethod.Method(0.43 ,     5, 0.01);

    //SimpleIteration 
    //SimpleIteration SimpleIter;
    //std::cout << "\nSolution Simple Iteration: " << SimpleIter.Iteration(1, 0.001);

    //Newton Method
   // NewtonMethod NewtonMeth;
    //std::cout <<"\nSolution Newton method for Polynominal: "<<NewtonMeth.Method(1, 2, 0.001);

    //Intersection Method
    //IntersectionMethod IntersectMeth;
    //std::cout <<"\nSolution Intersection method for Polynominal: "<<IntersectMeth.method(1, 2, 0.01);


}




