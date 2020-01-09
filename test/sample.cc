#include "fdlsgm.h"

#include<ctime>
#include<map>
#include<unordered_map>
#include<random>
#include<deque>

int
main(int argn, char** argv)
{
  std::unordered_map<size_t, fdlsgm::dls> pool;
  std::deque<fdlsgm::baseline> accumul, merged;

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
      const size_t n = n_elem*ig+ie;
      pool.emplace(n, a);
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
    const size_t n = n_elem*n_group + id;
    pool.emplace(n, a);
    printf("%3d ", -1); a.dprint();
  }
  printf("\n\n");

  /* first run */
  for (auto& v: pool) {
    bool flag = false;
    for (auto& bl: accumul) {
      const auto d = bl.argument(v.second)*180.0/M_PI;
      const auto l = bl.lateral_distance(v.second);
      const auto g = bl.gap_length(v.second);
      if (d < 45.0 && l < 3.0 && g < 0.5) {
        printf("# matched with [%08lx] (%6.4lf,%6.4f) / ", (size_t)&bl,d,l);
        flag = bl.append(v);
        if (flag) {
          printf("append to [%08lx].\n", (size_t)&bl);
        }
      } else {
        // printf("# not matched with [%08lx] (%6.4lf,%6.4f) \n",
        //        (size_t)&bl,d,l);
      }
    }
    if (accumul.empty() || flag == false) {
      printf("# create new baseline / ");
      fdlsgm::baseline baseline(v);
      accumul.push_back(baseline);
      printf("baseline [%08lx] pushed\n", (size_t)&accumul.back());
    }
  }
  /* iterate until convergence */
  bool changed = false;
  while (true) {
    for (auto& v: pool) {
      for (auto& bl: accumul) {
        const auto d = bl.argument(v.second)*180.0/M_PI;
        const auto l = bl.lateral_distance(v.second);
        const auto g = bl.gap_length(v.second);
        if (d < 45.0 && l < 3.0 && g < 0.5)
          changed |= bl.append(v);
      }
    }
    if (!changed) break;
    changed = false;
  }

  printf("# -------\n");

  while (!accumul.empty()) {
    const auto b = accumul.front();
    fdlsgm::baseline m(b);
    accumul.pop_front();
    const size_t N = accumul.size();
    for (size_t i=0; i<N; i++) {
      const auto c = accumul.front();
      accumul.pop_front();
      const auto d = b.argument(c)*180.0/M_PI;
      const auto l = b.lateral_distance(c);
      const auto g = b.gap_length(c);
      if ( d < 5.0 && l < 3.0 && g < 0.5) {
        auto k = m.size();
        m = merge_baseline(m,c);
        printf("merged!! (%ld = %ld + %ld)\n", m.size(), k, c.size());
      } else {
        accumul.push_back(c);
      }
    }
    if (m.size()>4)
      merged.push_back(m);
  }
  for (auto& e: merged) e.dprint();
  printf("# finally %ld baseline retrieved.\n", merged.size());

}
