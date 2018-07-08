#include "Matrices.h"

Matrix1D::Matrix1D(size_t size) : m_data(0., size)
{}

Matrix1D::Matrix1D(size_t size, fp value) : m_data(value, size)
{}

void Matrix1D::reset()
{ m_data.shift(m_data.size()); }

Matrix2D::Matrix2D(size_t width, size_t height) : Matrix2D(width, height, 0.)
{}

Matrix2D::Matrix2D(size_t width, size_t height, fp value) : width(width), height(height), m_data(value, width * height)
{
}

void Matrix2D::reset()
{ m_data.shift(m_data.size()); }

double& Matrix2D::operator()(int y, int x)
{ return m_data[y * width + x]; }

const double& Matrix2D::operator()(int y, int x) const
{ return m_data[y * width + x]; }

Matrix2D::data& Matrix2D::get()
{ return m_data; }

const Matrix2D::data& Matrix2D::get() const
{ return m_data; }

