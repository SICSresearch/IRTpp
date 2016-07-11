#ifndef GENERIC_HPP
#define GENERIC_HPP

#include <Types/Any.hpp>

namespace spgo
{

	template<typename T>
	class Generic: public Any
	{
	    public: 
	        Generic(T value) { data = new T(value); }
	        virtual ~Generic(){ delete (T*) data; }
	};
}

#endif