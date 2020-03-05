/**
 * @file sample_xyz.cc
 * @brief test in a 2000x2000x2000 space
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"

#include<ctime>
#include<map>
#include<unordered_map>
#include<random>
#include<deque>

int
main(int argn, char** argv)
{
  fdlsgm::accumulator<360> accumul;

  clock_t t = clock();
  std::mt19937 gen; gen.seed(t);
  std::uniform_real_distribution<double> location(-1000.0, 1000.0);
  std::normal_distribution<double> velocity(0.0, 2.0);
  std::normal_distribution<double> scatter(0.0, 0.1);

  const size_t n_elem       = 10;
  const size_t n_group      = 50;
  const size_t n_obstacle   = 50;
  const size_t n_distracter = 3000;

  clock_t begin = clock();
  for (size_t ig=0; ig<n_group; ig++) {
    const double x0 = location(gen);
    const double y0 = location(gen);
    const double z0 = location(gen);
    const double dx = velocity(gen);
    const double dy = velocity(gen);
    const double dz = velocity(gen);
    for (size_t ie=0; ie<n_elem; ie++) {
      const double x1 = x0+dx*ie+scatter(gen);
      const double y1 = y0+dy*ie+scatter(gen);
      const double z1 = z0+dz*ie+scatter(gen);
      const double x2 = x0+dx*(ie+1)+scatter(gen);
      const double y2 = y0+dy*(ie+1)+scatter(gen);
      const double z2 = z0+dz*(ie+1)+scatter(gen);
      accumul.insert(fdlsgm::dls({x1,y1,z1},{x2,y2,z2}));
    }
  }
  for (size_t ig=0; ig<n_obstacle; ig++) {
    const double x0 = location(gen);
    const double y0 = location(gen);
    const double z0 = location(gen);
    for (size_t ie=0; ie<n_elem; ie++) {
      const double x1 = x0+5*scatter(gen);
      const double y1 = y0+5*scatter(gen);
      const double z1 = z0+5*scatter(gen);
      const double x2 = x0+5*scatter(gen);
      const double y2 = y0+5*scatter(gen);
      const double z2 = z0+5*scatter(gen);
      accumul.insert(fdlsgm::dls({x1,y1,z1},{x2,y2,z2}));
    }
  }
  for (size_t id=0; id<n_distracter; id++) {
    const double x1 = location(gen);
    const double y1 = location(gen);
    const double z1 = location(gen);
    const double x2 = x1+velocity(gen);
    const double y2 = y1+velocity(gen);
    const double z2 = z1+velocity(gen);
    accumul.insert(fdlsgm::dls({x1,y1,z1},{x2,y2,z2}));
  }
  clock_t tic0 = clock();
  accumul.reallocate();
  clock_t tic1 = clock();
  accumul.merge();
  clock_t tic2 = clock();
  accumul.dprint(5);
  printf("## initialize: %lf sec\n", ((double)tic0-begin)/CLOCKS_PER_SEC);
  printf("## reallocate: %lf sec\n", ((double)tic1-tic0)/CLOCKS_PER_SEC);
  printf("## merge     : %lf sec\n", ((double)tic2-tic1)/CLOCKS_PER_SEC);
}
