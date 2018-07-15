#pragma once

#include <valarray>
#include "configStraylight.h"

class Matrix1D
{
public:
	using data = std::valarray<fp>;
	using reference = data::value_type&;
	using const_reference = const data::value_type&;
	using value_type = data::value_type;

	Matrix1D() = delete;

	Matrix1D(const Matrix1D& other) = delete;

	Matrix1D(Matrix1D&& other) = default;

	Matrix1D& operator=(const Matrix1D& other) = delete;

	Matrix1D& operator=(Matrix1D&& other) = default;

	explicit Matrix1D(size_t size);

	Matrix1D(size_t size, fp value);

	void reset();

	reference operator[](size_t pos)
	{ return m_data[pos]; }

	const_reference operator[](size_t pos) const
	{ return m_data[pos]; }

	data m_data;
};

class Matrix2D
{
public:
	using data = std::valarray<fp>;
	using reference = data::value_type&;
	using const_reference = const data::value_type&;
	using value_type = data::value_type;

	/*using iterator = decltype(data::begin());
	using const_iterator = data::const_iterator;*/
	Matrix2D() = delete;

	Matrix2D(const Matrix2D& other) = delete;

	Matrix2D(Matrix2D&& other) = default;

	Matrix2D& operator=(const Matrix2D& other) = delete;

	Matrix2D& operator=(Matrix2D&& other) = default;

	Matrix2D(size_t width, size_t height);

	Matrix2D(size_t width, size_t height, fp value);

	void clone(Matrix2D& rhs) const;

	reference operator()(int y, int x);

	const_reference operator()(int y, int x) const;

	data& get();

	const data& get() const;

	reference getLine(int y);

	reference operator[](size_t pos);

	const_reference operator[](size_t pos) const;

	void reset();
	/*iterator begin() noexcept { return m_data.begin(); }
	const_iterator begin() const noexcept { return m_data.begin(); }
	iterator end() noexcept { return m_data.end(); }
	const_iterator end() const noexcept { return m_data.end(); }*/
	// reference operator()(data::size_type x, data::size_type y) { return
	// m_data[pos]; }  const_reference operator()(data::size_type x,
	// data::size_type y) const { return m_data[pos]; }

	const size_t width;
	const size_t height;
	data m_data;
};