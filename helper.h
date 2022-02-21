#include <stdlib.h>
using namespace std;

template <typename T>
void info_w(T t) 
{
	std::cerr << t << std::endl ;
}

template<typename T, typename... Args>
void info_w(T t, Args... args) // recursive variadic function
{
	std::cerr << t << " " ;

	info_w(args...) ;
}

void info_w(){
	std::cerr << std::endl;
}

template<typename T, typename... Args>
void info(T t, Args... args) // recursive variadic function
{
	std::cerr << "INFO:\t"<< t ;

	info_w(args...) ;
}

template <typename T>
void error_w(T t) 
{
	std::cerr << t << std::endl ;
	exit(1);
}

template<typename T, typename... Args>
void error_w(T t, Args... args) // recursive variadic function
{
	std::cerr << t << " " ;

	error_w(args...) ;
}

void error_w(){
	std::cerr << endl;
	exit(1);
}

template<typename T, typename... Args>
void error(T t, Args... args) // recursive variadic function
{
	std::cerr << "ERROR:\t"<< t ;

	error_w(args...) ;
}

unsigned long long choose(unsigned long long n, unsigned long long k) {
	// stolen from Knuth
    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}