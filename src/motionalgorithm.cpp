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

using namespace SplineFittor;

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

int Motionalgorithm::getInflectionpointX_old(std::vector<double>pos_x, std::vector<double>pos_y, double &x_out)
{
	double maxPosionX = 0;
	double maxValue = 0;
	int n = pos_x.size();
	if (pos_x.size() != pos_y.size())
	{
		return -1;
	}
	std::vector<fitdata> data;
	data.resize(n);
	for (int i = 0; i < n; i++)
	{
		data[i].x = pos_x[i];
		data[i].y = pos_y[i];
	}
	std::sort(data.begin(), data.end(),
		[](const fitdata & a, const fitdata & b) {
		return a.x < b.x;
	});
	std::vector<double> x, y;
	x.resize(n);
	y.resize(n);
	for (int i = 0; i < n; ++i) {
		auto d = data.at(i);
		x[i] = d.x;
		y[i] = d.y;
	}
	std::vector<Coefficient> coe;
	Point p;
	p.xCoordinate = x;
	p.yCoordinate = y;
	p.f0 = 0;
	p.fn = 0;
	p.num = n;
	p.con = NOTAKNOT;
	spline(&p);
	coe = p.coe;
	int pointNum = 49;
	for (int i = 0; i < (n - 1); ++i) {
		double len = x[i + 1] - x[i];
		double fragment_len = len / (pointNum + 2 - 1);
		for (int j = 1; j <= pointNum; ++j) {
			double middlePoint_x = x[i] + j * fragment_len;
			double middlePoint_y = coe[i].a3 * j * fragment_len * j * fragment_len * j * fragment_len
				+ coe[i].b2 * j * fragment_len * j * fragment_len
				+ coe[i].c1 * j * fragment_len + coe[i].d0;
			maxValue < middlePoint_y ?
				(maxValue = middlePoint_y, maxPosionX = middlePoint_x) : 0;
		}
	}
	return maxPosionX;
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

int Motionalgorithm::spline(SplineFittor::Point *point)
{
	const std::vector<double> &x = point->xCoordinate;
	const std::vector<double> &y = point->yCoordinate;

	int n = point->num - 1;
	point->coe.resize(n);
	int i;
	Coefficient coe;
	Condition con = point->con;
	double *h, *d;
	double *a, *b, *c, *f, *m;
	double temp;

	if (n < 1)
	{
		return -2;
	}

	h = (double *)malloc((6 * n + 4) * sizeof(double));
	if (h == NULL)
	{
		return -1;
	}
	d = h + n;
	a = d + n;
	b = a + (n + 1);
	c = b + (n + 1);
	f = c + (n + 1);
	m = b;

	for (i = 0; i < n; i++)
	{
		h[i] = x[i + 1] - x[i];
		d[i] = (y[i + 1] - y[i]) / h[i];
	}
	a[0] = point->f0;
	c[n] = point->fn;
	switch (con)
	{
	case NATURAL:
		f[0] = 0;
		f[n] = 0;
		a[0] = 0;
		c[n] = 0;
		c[0] = 0;
		a[n] = 0;
		b[0] = 1;
		b[n] = 1;
		break;

	case CLAMPED:
		f[0] = 6 * (d[0] - a[0]);
		f[n] = 6 * (c[n] - d[n - 1]);
		a[0] = 0;
		a[n] = h[n - 1];
		c[n] = 0;
		c[0] = h[0];
		b[0] = 2 * h[0];
		b[n] = 2 * h[n - 1];
		break;

	case NOTAKNOT:
		f[0] = 0;
		f[n] = 0;
		a[0] = h[0];
		c[n] = h[n - 1];
		c[0] = -(h[0] + h[1]);
		a[n] = -(h[n - 2] + h[n - 1]);
		b[0] = h[1];
		b[n] = h[n - 2];
		break;
	}

	for (i = 1; i < n; i++)
	{
		a[i] = h[i - 1];
		b[i] = 2 * (h[i - 1] + h[i]);
		c[i] = h[i];
		f[i] = 6 * (d[i] - d[i - 1]);
	}

	if (n > 2)
	{
		for (i = 0; i <= n - 3; i++)
		{
			if (ABS(a[i + 1]) > ABS(b[i]))
			{
				SWAP(a[i + 1], b[i], temp);
				SWAP(b[i + 1], c[i], temp);
				SWAP(c[i + 1], a[i], temp);
				SWAP(f[i + 1], f[i], temp);
			}
			temp = a[i + 1] / b[i];
			a[i + 1] = 0;
			b[i + 1] = b[i + 1] - temp * c[i];
			c[i + 1] = c[i + 1] - temp * a[i];
			f[i + 1] = f[i + 1] - temp * f[i];
		}
	}
	if (n >= 2)
	{
		if (ABS(a[n - 1]) > ABS(b[n - 2]))
		{
			SWAP(a[n - 1], b[n - 2], temp);
			SWAP(b[n - 1], c[n - 2], temp);
			SWAP(c[n - 1], a[n - 2], temp);
			SWAP(f[n - 1], f[n - 2], temp);
		}
		if (ABS(c[n]) > ABS(b[n - 2]))
		{
			SWAP(c[n], b[n - 2], temp);
			SWAP(a[n], c[n - 2], temp);
			SWAP(b[n], a[n - 2], temp);
			SWAP(f[n], f[n - 2], temp);
		}
		temp = a[n - 1] / b[n - 2];
		a[n - 1] = 0;
		b[n - 1] = b[n - 1] - temp * c[n - 2];
		c[n - 1] = c[n - 1] - temp * a[n - 2];
		f[n - 1] = f[n - 1] - temp * f[n - 2];
		temp = c[n] / b[n - 2];
		c[n] = 0;
		a[n] = a[n] - temp * c[n - 2];
		b[n] = b[n] - temp * a[n - 2];
		f[n] = f[n] - temp * f[n - 2];
	}
	if (ABS(a[n]) > ABS(b[n - 1]))
	{
		SWAP(a[n], b[n - 1], temp);
		SWAP(b[n], c[n - 1], temp);
		SWAP(f[n], f[n - 1], temp);
	}
	temp = a[n] / b[n - 1];
	a[n] = 0;
	b[n] = b[n] - temp * c[n - 1];
	f[n] = f[n] - temp * f[n - 1];

	m[n] = f[n] / b[n];
	m[n - 1] = (f[n - 1] - c[n - 1] * m[n]) / b[n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		m[i] = (f[i] - (m[i + 2] * a[i] + m[i + 1] * c[i])) / b[i];
	}

	for (i = 0; i < n; i++)
	{
		coe.a3 = (m[i + 1] - m[i]) / (6 * h[i]);
		coe.b2 = m[i] / 2;
		coe.c1 = d[i] - (h[i] / 6) * (2 * m[i] + m[i + 1]);
		coe.d0 = y[i];
		point->coe[i] = coe;
	}
	free(h);
	return n + 1;
}