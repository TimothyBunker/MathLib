#include <stdio.h>
#include <stdlib.h>
#include <mathlib.h>
#include <constants.h>


double
expo(double base, int n) {
	int i;
	double total;
	total = 1.0;

	for(i = 0; i < n; ++i) {
		total *= base; 		
	}

	return total;
}

double
fact(int n) {
	if (n < 0){
		exit(1);
	}

	if (n == 0) { 
		return 1;
	}
	return n * fact(n - 1); 
}

double
sqrt_approx(double n){
	double epsilon, guess;

	if (n < 0) {
		return -1.0;
	}

	epsilon = 0.00001;
	guess = n / 2.0;

	while ((guess * guess) - n >= epsilon || (guess * guess) - n <= -epsilon) {
		guess = (guess + n / guess) / 2.0;
	}

	return guess;
}

double
rad_to_theta(double r) {
	double theta;
	theta = r * (180.0 / PI); 

	return theta;
}

double
theta_to_rad(double t) {
	double rad;
	rad = t * (PI / 180.0);

	return rad;
}

double
reduce_angle(double r) {
	while (r > 2 * PI) {
		r -= 2 * PI;
	} 
	while (r  < -2 * PI){
		r += 2 * PI;
	}

	return r;
}

double
mysin(double x, int rad){
	double sum;
	int terms, i;
	sum = 0.0;
	terms = 10;

	if (rad == 0) {
		x = theta_to_rad(x);
	}

	x = reduce_angle(x);
	for (i = 0; i < terms; ++i) {
		sum += (expo(-1.0, i) * expo(x, 2 * i + 1)) / fact(2 * i + 1);
	}

	return sum;

}

double
mycos(double x, int rad) {
	double sum;
	int terms, i;
	sum = 0.0;
	terms = 10;

	if (rad == 0) {
		x = theta_to_rad(x);
	}
	x = reduce_angle(x);
	for (i = 0; i < terms; ++i) {
		sum += (expo(-1.0, i) * expo(x, 2 * i)) / fact(2 * i); 
	}

	return sum;

}

double
mytan(double x, int rad) {
	return mysin(x, rad) / mycos(x, rad);
}

int
my_abs(int value) {
	return value < 0 ? -value : value;
}
