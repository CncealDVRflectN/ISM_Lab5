#pragma once
#include "Cloneable.h"

class PRNG : public Cloneable
{
public:
	virtual double next() const = 0;
	virtual double next(double from, double to) const = 0;
	virtual int nextInt(int from, int to) const = 0;
	virtual void reset() const = 0;
};

