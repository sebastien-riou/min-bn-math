//build with:
//gcc -std=c11 min-bn-math-test.c

#define MIN_BN_MATH_TEST 1
#include "min-bn-math.h"

int main(void){
	return min_bn_math_self_test();
}

