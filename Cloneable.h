#pragma once

class Cloneable
{
public:
	virtual Cloneable* clone() const = 0;
};