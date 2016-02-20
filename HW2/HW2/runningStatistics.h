#ifndef RUNNINGSTATS_H
#define RUNNINGSTATS_H

class runningStatistics
{
public:
	runningStatistics();
	void Clear();
	void Push(double x);
	long long NumDataValues() const;
	double Mean() const;
	double Variance() const;
	double StandardDeviation() const;
	double Skewness() const;
	double Kurtosis() const;

	friend runningStatistics operator+(const runningStatistics a, const runningStatistics b);
	runningStatistics& operator+=(const runningStatistics &rhs);

private:
	long long n;
	double M1, M2, M3, M4;
};

#endif
