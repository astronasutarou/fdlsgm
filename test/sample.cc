#include "fdlsgm.h"

#include<ctime>
#include<map>
#include<unordered_map>
#include<random>
#include<deque>

int
main(int argn, char** argv)
{
  fdlsgm::accumulator<80> accumul;

  clock_t t = clock();
  std::mt19937 gen; gen.seed(t);
  std::uniform_real_distribution<double> location(-1000.0, 1000.0);
  std::normal_distribution<double> velocity(0.0, 10.0);
  std::normal_distribution<double> scatter(0.0, 0.5);

  const size_t n_elem       = 10;
  const size_t n_group      = 50;
  const size_t n_obstacle   = 50;
  const size_t n_distracter = 3000;

  clock_t begin = clock();
  for (size_t ig=0; ig<n_group; ig++) {
    const double x0 = location(gen);
    const double y0 = location(gen);
    const double dx = velocity(gen);
    const double dy = velocity(gen);
    for (size_t ie=0; ie<n_elem; ie++) {
      const double x = x0+dx*ie, y = y0+dy*ie, z = (double)ie;
      const double x1 = scatter(gen), y1 = scatter(gen);
      const double x2 = scatter(gen), y2 = scatter(gen);
      accumul.insert(fdlsgm::dls({x1+x,y1+y,z},{x2+x+dx,y2+y+dy,z+1}));
    }
  }
  for (size_t ig=0; ig<n_obstacle; ig++) {
    const double x = location(gen);
    const double y = location(gen);
    for (size_t ie=0; ie<n_elem; ie++) {
      const double z = (double)ie;
      const double x1 = scatter(gen), y1 = scatter(gen);
      const double x2 = scatter(gen), y2 = scatter(gen);
      accumul.insert(fdlsgm::dls({x1+x,y1+y,z},{x2+x,y2+y,z+1}));
    }
  }
  for (size_t id=0; id<n_distracter; id++) {
    const double x0 = location(gen);
    const double y0 = location(gen);
    const double z0 = (location(gen)+1000.)/2000.*n_elem;
    const double dx = velocity(gen);
    const double dy = velocity(gen);
    accumul.insert(fdlsgm::dls({x0,y0,z0},{x0+dx,y0+dy,z0+1.0}));
  }
  clock_t tic0 = clock();
  accumul.reallocate();
  clock_t tic1 = clock();
  accumul.coalesce();
  clock_t tic2 = clock();
  accumul.dprint(5);
  printf("## initialize: %lf sec\n", ((double)tic0-begin)/CLOCKS_PER_SEC);
  printf("## reallocate: %lf sec\n", ((double)tic1-tic0)/CLOCKS_PER_SEC);
  printf("## coalesce  : %lf sec\n", ((double)tic2-tic1)/CLOCKS_PER_SEC);
}
