#include "stdafx.h"

#include <iostream>
#include <fstream>


using namespace std;

int conv(double *, int, double *, int, double *);
double f(double);
double h(double);
void initialCond(double *);
void mattimes(double t, double a[][3], double m[][3]);
void matsub(double a[][3], double b[][3], double c[][3]);
void matvecmult(double m[][3], double *v, double *prod);


int main()
{
	const int SIZE = 10001;
	double arrf[SIZE];
	double arrh[SIZE];
	double arrfh[SIZE*2-1];
	double arry0[SIZE];
	double arry[SIZE*2-1];
	ofstream outfile("Problem3_2");
	ofstream outfilef("f");
	ofstream outfileh("h");
	ofstream outfilefh("fh");

	//create f(t)
	for (int i = 0; i < SIZE; i ++) {
		arrf[i] = f(i*0.001);
	}

	//create h(t)
	for (int i = 0; i < SIZE; i++) {
		arrh[i] = h(i*0.001);
		//outfile << (t*0.001) << " " << arry[t] << endl;
	}

	//convolute f and h to create fh
	conv(arrf, SIZE, arrh, SIZE, arrfh);
	
	//create y0
	initialCond(arry0);
	//for (int i = 0; i < SIZE; i++) {
		//outfile << (i*0.001) << " " << arry0[i] << endl;
	//}

	for (int t = 0; t < SIZE; t++) {
		arry[t] = arrfh[t]; + arry0[t];
		outfile << (t*0.001) << "   " << arry[t] << endl;
		//outfile << t*0.001 << " " << (arry[t]);
		cout << arry[t] << endl;
	}
	for (int t = 0; t < SIZE; t++) {
		arry[t+SIZE] = arrfh[t+SIZE];
		//outfile << (t*0.001) << "   " << arry[t] << endl;
		outfile << t*0.001+10.001 << " " << (arry[t+SIZE])<<endl;
		cout << arry[t] << endl;
	}

    return 0;
}

int conv(double *f1, int len1, double *f2, int len2, double* y) {
	int totLen = len1 + len2 - 1;
	int val, val2;
	double temp;
	for (int n = 0; n < totLen; n++) {
		temp = 0;
		if (n >= len2) { val = n - len2 + 1; }
		else { val = 0; }
		if (n >= len1) { val2 = len1 - 1; }
		else { val2 = n; }
		for (int i = val; i <= val2; i++) {
			temp += (f1[i] * f2[n - i]);
		}
		y[n] = temp;
	}
	return totLen;
}
double f(double t) {
	return sin(4 * 3.14159 * t);
}

double h(double t) {
	return (-2 / 153)*exp(-2 * t) + 0.16392*exp(-0.5*t)*cos(6 * t - 85.4261);
}

void initialCond(double * arry0){
double m1[3][3];
double m2[3][3];
double I[3][3] = { { 1,0,0 },{ 0,1,0 },{ 0,0,1 } };

double deltat, a1, a2, a3, t1, t2, x[3];
int n;

cout << "Enter initial condition y0(0) = ";
cin >> x[0];
cout << "Enter initial condition y1(0) = ";
cin >> x[1];
cout << "Enter initial condition y2(0) = ";
cin >> x[2];
cout << "Enter number of iterations: ";
cin >> n;
cout << "Enter min bound of t: ";
cin >> t1;
cout << "Enter max bound of t: ";
cin >> t2;
cout << "For x'3(t) + a2x3(t) + a1x2(t) + a0x1(t) = 0" << endl;
cout << "Enter value for a0 in (D-a)y0(t)=0 : ";
cin >> a1;
cout << "Enter value for a1 in (D-a)y0(t)=0 : ";
cin >> a2;
cout << "Enter value for a2 in (D-a)y0(t)=0 : ";
cin >> a3;

double a[3][3] = { { 0,-1,0 },{ 0,0,-1 },{ a1,a2,a3 } };
deltat = (t2 - t1) / n;



//outfile << "0" << " " << x[0] << endl;
	arry0[0] = x[0];
for (int i = 1; i < n; i++) {

	mattimes(deltat, a, m1);
	matsub(I, m1, m2);


	// x = (1 + m1)*x;
	matvecmult(m2, x, x);
	arry0[i] = x[0];
}






}



void mattimes(double t, double a[][3], double m[][3])
{
	int i, j;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			m[i][j] = t*a[i][j];
}
void matsub(double a[][3], double b[][3], double c[][3])
{
	int i, j;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			c[i][j] = a[i][j] - b[i][j];
}
void matvecmult(double m[][3], double *v, double *prod)
{
	double sum;
	int i, j;
	for (i = 0; i < 3; i++) {
		sum = 0;
		for (j = 0; j < 3; j++) {
			sum += m[i][j] * v[j];
		}
		prod[i] = sum;
	}
}

//#include<iostream>
//
//using namespace std;
//
//int conv(double *, int, double *, int, double *);
//
//int main() {
//	double leny;
//	double y[100] = { 0 };
//	double f1[] = { 0,1,2,3,2,1 };
//	int len1 = 6;
//	double f2[] = { -2,-2,-2,-2,-2,-2,-2 };
//	int len2 = 7;
//	double f3[] = { 1,-1,1,-1 };
//	int len3 = 4;
//	double f4[] = { 0,0,0,-3,-3 };
//	int len4 = 5;
//
//	leny = conv(f1, len1, f1, len1, y);
//	cout << "f1*f1 = ";
//	for (int i = 0; i < leny; i++) {
//		cout << y[i] << " ";
//	}
//	cout << endl;
//	leny = conv(f1, len1, f2, len2, y);
//	cout << "f1*f2 = ";
//	for (int i = 0; i < leny; i++) {
//		cout << y[i] << " ";
//	}
//	cout << endl;
//	leny = conv(f1, len1, f3, len3, y);
//	cout << "f1*f3 = ";
//	for (int i = 0; i < leny; i++) {
//		cout << y[i] << " ";
//	}
//	cout << endl;
//	leny = conv(f2, len2, f3, len3, y);
//	cout << "f2*f3 = ";
//	for (int i = 0; i < leny; i++) {
//		cout << y[i] << " ";
//	}
//	cout << endl;
//	leny = conv(f1, len1, f4, len4, y);
//	cout << "f1*f4 = ";
//	for (int i = 0; i < leny; i++) {
//		cout << y[i] << " ";
//	}
//	cout << endl;
//}