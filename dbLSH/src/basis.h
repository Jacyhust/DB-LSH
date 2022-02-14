/**
 * @file basis.h
 *
 * @brief A set of basic tools.
 */
#pragma once
#include <string>
#include <iostream>
#include <time.h>
#include <vector>
namespace lsh
{
	class progress_display
	{
	public:
		explicit progress_display(
			unsigned long expected_count,
			std::ostream& os = std::cout,
			const std::string& s1 = "\n",
			const std::string& s2 = "",
			const std::string& s3 = "")
			: m_os(os), m_s1(s1), m_s2(s2), m_s3(s3)
		{
			restart(expected_count);
		}
		void restart(unsigned long expected_count)
		{
			_count = _next_tic_count = _tic = 0;
			_expected_count = expected_count;
			m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
				<< m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
				<< std::endl
				<< m_s3;
			if (!_expected_count)
			{
				_expected_count = 1;
			}
		}
		unsigned long operator += (unsigned long increment)
		{
			if ((_count += increment) >= _next_tic_count)
			{
				display_tic();
			}
			return _count;
		}
		unsigned long  operator ++ ()
		{
			return operator += (1);
		}
		unsigned long count() const
		{
			return _count;
		}
		unsigned long expected_count() const
		{
			return _expected_count;
		}
	private:
		std::ostream& m_os;
		const std::string m_s1;
		const std::string m_s2;
		const std::string m_s3;
		unsigned long _count, _expected_count, _next_tic_count;
		unsigned _tic;
		void display_tic()
		{
			unsigned tics_needed = unsigned((double(_count) / _expected_count) * 50.0);
			do
			{
				m_os << '*' << std::flush;
			} while (++_tic < tics_needed);
			_next_tic_count = unsigned((_tic / 50.0) * _expected_count);
			if (_count == _expected_count)
			{
				if (_tic < 51) m_os << '*';
				m_os << std::endl;
			}
		}
	};
	/**
	 * A timer object measures elapsed time, and it is very similar to boost::timer.
	 */
	class timer
	{
	public:
		timer() : time(double(clock())) {};
		~timer() {};
		/**
		 * Restart the timer.
		 */
		void restart()
		{
			time = double(clock());
		}
		/**
		 * Measures elapsed time.
		 *
		 * @return The elapsed time
		 */
		double elapsed()
		{
			return (double(clock()) - time) / CLOCKS_PER_SEC;
		}
	private:
		double time;
	};
}

float cal_inner_product(float* v1, float* v2, int dim);

float cal_dist(float* v1, float* v2, int dim);

float cal_l2_dist(std::vector<float>& v1, std::vector<float>& v2);

float cal_l2_dist(const std::vector<float>& v1, const std::vector<float>& v2);

void set_rmin(std::string& datatsetName, float& R_min);

std::vector<std::vector<float>> Random_pivot(float** raw_data_, int pivot_num_, int N_, int dim, int repeat_num);

float distance_Between_Piovt_Vector(std::vector<std::vector<float>>& pivot_vector_);

template <class T>
void clear_2d_array(T** array, int n)
{
	for (int i = 0; i < n; ++i) {
		delete[] array[i];
	}
	delete[] array;
}

void showMemoryInfo();

//void lshknn();

//bool Is_Intersect(float*& mbr, float*& data, int& dim);