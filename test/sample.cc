#include "fdlsgm.h"

#include<ctime>
#include<map>
#include<unordered_map>
#include<random>
#include<deque>

int
main(int argn, char** argv)
{
  fdlsgm::accumulator<45> accumul;

  clock_t t = clock();
  std::mt19937 gen; gen.seed(t);
  std::uniform_real_distribution<double> location(-100.0, 100.0);
  std::normal_distribution<double> velocity(0.0, 1.0);
  std::normal_distribution<double> scatter(0.0, 0.3);

  const size_t n_elem     = 10;
  const size_t n_group    = 50;
  const size_t n_distract = 300;

  for (size_t ig=0; ig<n_group; ig++) {
    const double x0 = location(gen);
    const double y0 = location(gen);
    const double dx = velocity(gen);
    const double dy = velocity(gen);
    for (size_t ie=0; ie<n_elem; ie++) {
      const double x = x0+dx*ie, y = y0+dy*ie, z = (double)ie;
      const double x1 = scatter(gen), y1 = scatter(gen);
      const double x2 = scatter(gen), y2 = scatter(gen);
      fdlsgm::dls
        a({x1+x,y1+y,z},{x2+x+dx,y2+y+dy,z+1});
      accumul.insert(a);
      printf("%3ld ",1+(ig%8)); a.dprint();
    }
  }

  for (size_t id=0; id<n_distract; id++) {
    const double x0 = location(gen);
    const double y0 = location(gen);
    const double z0 = (location(gen)+100.)/200.*n_elem;
    const double dx = velocity(gen);
    const double dy = velocity(gen);
    fdlsgm::dls a({x0,y0,z0},{x0+dx,y0+dy,z0+1.0});
    accumul.insert(a);
    printf("%3d ", -1); a.dprint();
  }
  printf("\n\n");

  printf("N(segment) : %ld\n", accumul.count_dls());
  printf("N(baseline): %ld\n", accumul.count_baseline());
}
