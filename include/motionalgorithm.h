#pragma once

#define DLL_EXPORTS extern "C" __declspec(dllexport)
#define DLL_CLASSEXP __declspec(dllexport)
#include "dataStructure.h"
#include <vector>

namespace SplineFittor {

	#define ABS(p) ( ((p) > 0) ? (p) : -(p) )
	#define SWAP(x, y, temp) (temp) = (y); (y) = (x); (x) = (temp);

	struct fitdata
	{
		double x;
		double y;
	};
	struct Coefficient
	{
		double a3;
		double b2;
		double c1;
		double d0;
	};
	enum Condition
	{
		NATURAL,
		CLAMPED,
		NOTAKNOT
	};
	struct Point
	{
		std::vector<Coefficient> coe;
		std::vector<double> xCoordinate;
		std::vector<double> yCoordinate;
		double f0;
		double fn;
		int num;
		Condition con;
	};
}

class DLL_CLASSEXP Motionalgorithm
{
public:
	Motionalgorithm();
	~Motionalgorithm();

	//最小二乘法拟合法寻找拐点坐标(可用于计算对焦清晰位置)
	int getInflectionpointX(std::vector<double>pos_x,std::vector<double>pos_y,double &x_out);

	//最小二乘法拟合法寻找拐点坐标(可用于计算对焦清晰位置)
	int getInflectionpointX_old(std::vector<double>pos_x, std::vector<double>pos_y, double &x_out);

private:
	double vertex_x(const double a, const double b);
	double parabola_value(const double a, const double b, const double c, double x);
	int polyfit(const double *x, const double *y, int xyLength, int poly_n, std::vector<double> &out_factor, double &out_chisq);
	int spline(SplineFittor::Point *point);
};

