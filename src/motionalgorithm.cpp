#include "commitid.h"
#include "alg_base_common.h"
#include "motionalgorithm.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

Motionalgorithm::Motionalgorithm()
{

}

Motionalgorithm::~Motionalgorithm()
{

}

int Motionalgorithm::getInflectionpointX(std::vector<double>pos_x, std::vector<double>pos_y, double &x_out)
{
	struct fitdata
	{
		double x;
		double y;
	};
	int size = pos_x.size();
	if (pos_x.size()!= pos_y.size())
	{
		return -1;
	}
	std::vector<fitdata> datas;
	datas.resize(size);
	for (int n=0;n<size;n++)
	{
		datas[n].x = pos_x[n];
		datas[n].y = pos_y[n];
	}
	std::sort(datas.begin(), datas.end(),
		[](const fitdata & a, const fitdata & b) {
		return a.x < b.x;
	});

	int maxIndex = 0;
	for (int i = 1; i < datas.size(); ++i) {
		if (datas[i].y > datas[maxIndex].y) {
			maxIndex = i;
		}
	}

	double x[5];
	double y[5];
	if (datas.size() > 10) {
		for (int i = 0; i < 5; i++) {
			if ((maxIndex + i - 2) < 0 || (maxIndex + i - 2) > (datas.size() - 1)) {
				return datas[maxIndex].x;
			}
			y[i] = datas[maxIndex + i - 2].y;
			x[i] = datas[maxIndex + i - 2].x;
		}
	}
	else {
		return datas[maxIndex].x;
	}

	int len = 5;
	int poly_n = 2;
	std::vector<double> out_factor;
	double out_chisq;
	polyfit(x, y, len, poly_n, out_factor, out_chisq);//out_factor 2-a, 1-b, 0-c

	double a = out_factor[2];
	double b = out_factor[1];
	double c = out_factor[0];

	double vertex_x_value = vertex_x(a, b);
	double vertex_y_value = parabola_value(a, b, c, vertex_x_value);
	return vertex_x_value;

	return 0;
}

int Motionalgorithm::polyfit(const double *x, const double *y, int xyLength, int poly_n, std::vector<double> &out_factor, double &out_chisq)
{
	/*
	 * x：自变量，视差
	 * y：因变量，距离
	 * xyLength: x、y长度
	 * poly_n：拟合的阶次
	 * out_factor：拟合的系数结果，从0阶到poly_n阶的系数
	 * out_chisq：拟合曲线与数据点的优值函数最小值 ,χ2 检验
	*/

	gsl_matrix *XX = gsl_matrix_alloc(xyLength, poly_n + 1);
	gsl_vector *c = gsl_vector_alloc(poly_n + 1);
	gsl_matrix *cov = gsl_matrix_alloc(poly_n + 1, poly_n + 1);
	gsl_vector *vY = gsl_vector_alloc(xyLength);

	for (size_t i = 0; i < xyLength; i++) {
		gsl_matrix_set(XX, i, 0, 1.0);
		gsl_vector_set(vY, i, y[i]);
		for (int j = 1; j <= poly_n; j++) {
			gsl_matrix_set(XX, i, j, pow(x[i], j));
		}
	}
	gsl_multifit_linear_workspace *workspace = gsl_multifit_linear_alloc(xyLength, poly_n + 1);
	int r = gsl_multifit_linear(XX, vY, c, cov, &out_chisq, workspace);
	gsl_multifit_linear_free(workspace);
	out_factor.resize(c->size, 0);
	for (size_t i = 0; i < c->size; ++i) {
		out_factor[i] = gsl_vector_get(c, i);
	}

	gsl_vector_free(vY);
	gsl_matrix_free(XX);
	gsl_matrix_free(cov);
	gsl_vector_free(c);

	return r;
}

double Motionalgorithm::vertex_x(const double a, const double b)
{
	return -b / (2.0 * a);
}

double Motionalgorithm::parabola_value(const double a, const double b, const double c, double x)
{
	return a * x * x + b * x + c;
}