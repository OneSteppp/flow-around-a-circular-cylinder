#include<iostream>
#include<fstream>
#include<cmath>
#include <iomanip>
using namespace std;

double* matrix_subtract(double *a, double *b,int &m)
{
	for (int i = 0; i < m; i++)
	{
		b[i] = a[i] - b[i];
	}
	return b;
}//矩阵相减
double Norm_2(double *a, int &m)
{
	double n = 0;
	for (int i = 0; i < m; i++)
	{
		n = n + a[i] * a[i];
	}
	n = sqrt(n);
	return n;
}//求范数

void Gauss_Seidel(double *A, double *y, double *r, int &m, double error)
{
	double *y_old;
	y_old = new double[m];
	double N_2;
	do
	{
		for (int i = 0; i < m; i++)
		{
			y_old[i] = y[i];
		}
		for (int i = 0; i < m; i++)
		{
			double sum = 0;
			for (int j = 0; j < m; j++)
			{
				if (j == i) continue;//跳过Aii
				sum = sum + A[i*m + j] * y[j];
			}
			y[i] = (r[i] - sum) / A[i*m + i];
		}
		N_2 = Norm_2(matrix_subtract(y, y_old, m), m);
		//cout << N_2<<'\t';
	} while (N_2> error);
	delete[] y_old;
	y_old = NULL;
}//高斯-赛德尔迭代法

int main()
{
	int m, n;//横纵网格数
	cout << "请输入横竖坐标网格数：";
	cin >> m >> n;
	double delt_Xi = 4.0 / m;
	double delt_phi = 2.0 / n;//网格的长度和宽度
	double delt_t;
	cout << "请输入时间步长:";
	cin >> delt_t;
	double *Xi, *phi;//横坐标与纵坐标
	double *y, *Y_old, *A, *r;//A为系数矩阵
	int size = (m + 1)*(n + 1);//系数矩阵大小
	double y_Xi, y_phi, y_Xi_phi, y_Xi_Xi, y_phi_phi;//偏导数

	Xi = new double[m+1];
	phi = new double[n+1];
	y = new double[size];
	Y_old = new double[size];
	A = new double[size*size];
	r = new double[size];
	for (int j = 0; j <= m; j++)
	{
		Xi[j] = -4 + j*delt_Xi;
	}//横坐标赋值[-4:0.1:0]
	for (int i = 0; i <= n; i++)
	{
		phi[i] = i*delt_phi;
	}//纵坐标赋值[0:0.1:2]
	//
	//
	for (int i = 0; i < size; i++)
	{
		y[i] = 0;
		r[i] = 0;
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			A[i*size + j] = 0;
		}
	}//矩阵赋初值
	//
	//
	for (int j = 0.75*m; j <= m; j++)
	{
		int i = 0;
		y[i*(m + 1) + j] = sqrt(1 - Xi[j] * Xi[j]);
	}//下边界
	for (int j = 0; j <= m; j++)
	{
		int i = n;
		y[i*(m + 1) + j] = 2;
	}//上边界
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			y[i*(m + 1) + j] = y[j]+(y[n*(m+1)+j]-y[j])*i/n;
		}
	}//插值得到场中间值
	/*ofstream zerofile("y_0.txt");
	for (int j = 0; j <= m; j++)
	{
		for (int i = 0; i <= n; i++)
		{
			zerofile << setprecision(5) << y[i*(m + 1) + j] << '\t';
		}
		cout << '\n';
	}
	zerofile.close();*/
	//
	//
	for (int j = 0; j < 0.75*m; j++)
	{
		int i = 0;
		A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1;
		r[j] = 0;
	}
	for (int j = 0.75*m; j <= m; j++)
	{
		int i = 0;
		A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1;
		r[i*(m + 1) + j] = sqrt(1 - Xi[j] * Xi[j]);
	}//y[i,j] i=0 j [0 40] 下边界
	for (int j = 0; j <= m; j++)
	{
		int i = n;
		A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1;
		r[i*(m + 1) + j] = 2;
	}//y[i,j] i=20 j [0 40]  上边界
	for (int i = 0; i <= n; i++)
	{
		int j = 0;
		//A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1;
		//A[(i*(m + 1) + j)*size + i*(m + 1) + j + 1] = -1;
		//r[i*(m + 1) + j] = 0;
		A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1;
		r[i*(m + 1) + j] = phi[i];
	}//y[i,j] i [0 20] j=0  左边界
	for (int i = 0; i <= n; i++)
	{
		int j = m;
		A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1;
		A[(i*(m + 1) + j)*size + i*(m + 1) + j - 1] = -1;
		r[i*(m + 1) + j] = 0;
	}//y[i,j] i [1 19] j=40 右边界


	double norm_2;//2-范数
	do
	{
		for (int i = 0; i < size; i++)
		{
			Y_old[i] = y[i];
		}
		//将y赋值给Y_old

		for (int i = 1; i <= n-1; i++)
		{
			for (int j = 1; j <= m-1; j++)
			{
				y_phi = 0.5*(y[(i + 1)*(m + 1) + j] - y[(i - 1)*(m + 1) + j]) / delt_phi;
				y_Xi = 0.5*(y[i *(m + 1) + j + 1] - y[i*(m + 1) + j - 1]) / delt_Xi;
				y_Xi_phi = 0.25*(y[(i + 1) *(m + 1) + j + 1] - y[(i + 1) *(m + 1) + j - 1] - y[(i - 1) *(m + 1) + j + 1] + y[(i - 1) *(m + 1) + j - 1]) / (delt_phi*delt_Xi);
				y_phi_phi = (y[(i + 1)*(m + 1) + j] - 2 * y[i *(m + 1) + j] + y[(i - 1)*(m + 1) + j]) / (delt_phi*delt_phi);
				//
				A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1 / delt_t + 2 * y_phi*y_phi / (delt_Xi*delt_Xi);
				A[(i*(m + 1) + j)*size + i*(m + 1) + j - 1] = -y_phi*y_phi / (delt_Xi*delt_Xi);
				A[(i*(m + 1) + j)*size + i*(m + 1) + j + 1] = -y_phi*y_phi / (delt_Xi*delt_Xi);
				//
				r[i*(m + 1) + j] = y[i*(m + 1) + j] / delt_t - y_phi*y_Xi*y_Xi_phi;
				//r[i*(m + 1) + j] = y[i*(m + 1) + j] / delt_t - 2 * y_phi*y_Xi*y_Xi_phi + (1 + y_Xi*y_Xi)*y_phi_phi;
				//cout << A[(i*(m + 1) + j)*size + i*(m + 1) + j] << '\t' << A[(i*(m + 1) + j)*size + i*(m + 1) + j - 1] << '\t' << A[(i*(m + 1) + j)*size + i*(m + 1) + j + 1] << '\t' << r[i*(m + 1) + j]<<endl;
			}
		}
		//
		Gauss_Seidel(A, y, r, size, 1e-6);//解A*y=r；
		for (int i = 1; i <= n - 1; i++)
		{
			for (int j = 1; j <= m - 1; j++)
			{
				A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 0;
				A[(i*(m + 1) + j)*size + i*(m + 1) + j - 1] = 0;
				A[(i*(m + 1) + j)*size + i*(m + 1) + j + 1] = 0;
				r[i*(m + 1) + j] = 0;
			}
		}//清除A矩阵内容
		/*ofstream onefile("y_1.txt");
		for (int j = 0; j <= m; j++)
		{
			for (int i = 0; i <= n; i++)
			{
				onefile << setprecision(5) << y[i*(m + 1) + j] << '\t';
			}
			cout << '\n';
		}
		onefile.close();
		//cout << y[345];
		cout << endl;*/
		//
		for (int i = 1; i <= n-1; i++)
		{
			for (int j = 1; j <= m-1; j++)
			{
				y_phi = 0.5*(y[(i + 1)*(m + 1) + j] - y[(i - 1)*(m + 1) + j]) / delt_phi;
				y_Xi = 0.5*(y[i *(m + 1) + j + 1] - y[i*(m + 1) + j - 1]) / delt_Xi;
				y_Xi_phi = 0.25*(y[(i + 1) *(m + 1) + j + 1] - y[(i + 1) *(m + 1) + j - 1] - y[(i - 1) *(m + 1) + j + 1] + y[(i - 1) *(m + 1) + j - 1]) / (delt_phi*delt_Xi);
				y_Xi_Xi = (y[i *(m + 1) + j + 1] - y[i *(m + 1) + j] + y[i *(m + 1) + j - 1]) / (delt_Xi*delt_Xi);
				//
				A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 1 / delt_t + 2 * (1 + y_Xi*y_Xi) / (delt_phi*delt_phi);
				A[(i*(m + 1) + j)*size + (i - 1)*(m + 1) + j] = -(1 + y_Xi*y_Xi) / (delt_phi*delt_phi);
				A[(i*(m + 1) + j)*size + (i + 1)*(m + 1) + j] = -(1 + y_Xi*y_Xi) / (delt_phi*delt_phi);
				//
				r[i*(m + 1) + j] = y[i*(m + 1) + j] / delt_t - y_phi*y_Xi*y_Xi_phi;
				//r[i*(m + 1) + j] = y[i*(m + 1) + j] / delt_t - 2 * y_phi*y_Xi*y_Xi_phi + y_phi*y_phi*y_Xi_Xi;
				//cout << A[(i*(m + 1) + j)*size + i*(m + 1) + j] << '\t' << A[(i *(m + 1) + j)*size + (i + 1)*(m + 1) + j] << '\t' << A[(i *(m + 1) + j)*size + (i - 1)*(m + 1) + j] << '\t' << r[i*(m + 1) + j] << endl;
			}
			//Gauss_Seidel(A, y, r, size, 1e-6);//解A*y=r；//解A*y=r;
		}
		//
		Gauss_Seidel(A, y, r, size, 1e-6);//解A*y=r；//解A*y=r;
		for (int i = 1; i <= n - 1; i++)
		{
			for (int j = 1; j <= m - 1; j++)
			{
				A[(i*(m + 1) + j)*size + i*(m + 1) + j] = 0;
				A[(i*(m + 1) + j)*size + (i - 1)*(m + 1) + j] = 0;
				A[(i*(m + 1) + j)*size + (i + 1)*(m + 1) + j] = 0;
				r[i*(m + 1) + j] = 0;
			}
		}//清除A矩阵内容

		/*ofstream twofile("y_2.txt");
		for (int j = 0; j <= m; j++)
		{
			for (int i = 0; i <= n; i++)
			{
				twofile << setprecision(5) << y[i*(m + 1) + j] << '\t';
			}
			cout << '\n';
		}
		twofile.close();
		cout << endl;*/
		//
		norm_2 = Norm_2(matrix_subtract(y, Y_old, size), size);
		/*ofstream threefile("y_3.txt");
		for (int j = 0; j <= m; j++)
		{
			for (int i = 0; i <= n; i++)
			{
				threefile << setprecision(5) << Y_old[i*(m + 1) + j] << '\t';
			}
			cout << '\n';
		}
		threefile.close();
		cout << endl;*/
		cout << norm_2 << endl;
		//break;
	} while (norm_2 > 1e-6);
	ofstream Ofile("y_result.txt");
	for (int j = 0; j <= m ; j++)
	{
		Ofile << Xi[j] << '\t';
		for (int i = 0; i <= 5; i++)
		{
			Ofile  << setprecision(5) << y[n*i*(m + 1)/5 + j] << '\t';
		}
		cout << '\n';
	}
	Ofile.close();

	delete[] A;
	A = NULL;
	delete[] Y_old;
	Y_old = NULL;
	delete[] r;
	r = NULL;
	delete[] y;
	y = NULL;
	delete[] Xi;
	Xi = NULL;
	delete[] phi;
	phi = NULL;
	return 0;
}