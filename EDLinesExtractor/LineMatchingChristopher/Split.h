#pragma once
#include <cstddef>



struct split1
{
	enum empties_t { empties_ok, no_empties };
};

template <typename Container>
Container& split(
	Container&                            result,
	const typename Container::value_type& s,
	const typename Container::value_type& delimiters,
	split1::empties_t                      empties = split1::empties_ok )
{
	result.clear();
	size_t current;
	size_t next = -1;
	do
	{
		if (empties == split1::no_empties)
		{
			next = s.find_first_not_of( delimiters, next + 1 );
			if (next == Container::value_type::npos) break;
			next -= 1;
		}
		current = next + 1;
		next = s.find_first_of( delimiters, current );
		result.push_back( s.substr( current, next - current ) );
	}
	while (next != Container::value_type::npos);
	return result;
}
