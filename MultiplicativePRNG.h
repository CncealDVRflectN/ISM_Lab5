#pragma once
#include "PRNG.h"

class MultiplicativePRNG : public PRNG
{
private:
	const long long module;
	const long long seed;
	const int multiplier;

	mutable long long last;
public:
	MultiplicativePRNG(long long module, long long seed, int multiplier);
	MultiplicativePRNG(const MultiplicativePRNG* source);
	~MultiplicativePRNG();

	double next() const override;
	double next(double from, double to) const override;
	int nextInt(int from, int to) const override;
	void reset() const override;
	MultiplicativePRNG* clone() const override;
};

