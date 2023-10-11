#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double generate_uniform_random(void)
/* Generate a uniform random number between 0 and 1, inclusive */
/* Note: this is NOT the best random number out there, but it's quick and simple */
{
{
  double dzeta;
  //srand(time(0));
  double Bottom = -sqrt(3);
  dzeta = Bottom + (((double)rand() / (double)RAND_MAX)*2*sqrt(3));
  return dzeta;
}
}

double gauss(void)
/* Generate a pseudo-random number from a near-Gaussian distribution with zero average and unit variance */
/* Central Limit Theorem: sum of 12 uniform random numbers [0,1] has average of 6 and variance of 1 */
{
  double sum = -6.0;
  for (int i = 0; i < 12; i++)
  {
    sum += generate_uniform_random();
  }
  return sum;
}
