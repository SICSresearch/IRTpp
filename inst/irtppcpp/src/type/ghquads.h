#ifndef GHQUADS_H_
#define GHQUADS_H_

namespace irtpp
{

  extern const double full_gh_quads[];
  extern const double full_gh_weigths[];
  const double * quads(int i);
  const double * weights(int i);

}

#endif
