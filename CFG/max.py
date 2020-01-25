from sympy import Max, Min, symbols
from functools import reduce, lru_cache
from scipy.stats.distributions import binom

def mmax(xs):
  return reduce(Max, xs)

def compositions(n, k):
  if n < 0 or k < 0: return
  if n == 0 and k == 0:
    yield ()
    return

  for j in range(0, n+1):
    for xs in compositions(n-j, k-1):
      yield (j,) + xs

# Lower bound.
def OPT(n, p):
  X = binom(n, p)
  return sum(
    X.pmf(i) * sum(
      2 * Min(e, o) + 2
      for e, o in compositions(i, 2) if 2*e <= i and 2*o <= i+1
    ) for i in range(0, n+1)
  )


@lru_cache(maxsize=None)
def V(state, n, p):
  if n % 2 == 1 or n < 4: return 0
  elif state == 'L': return 1*p**2 + V('L',n-2,p)
  elif state == 'R':
    return mmax([
      V('S',a,p) + V('S',b,p)
      for a,b in compositions(n,2)
    ])
  elif state == 'S':
    return mmax([
      1*p**2 + V('S',n-2,p),
      3*p**2 + mmax([V('L',a,p) + V('S',b,p) for a,b in compositions(n-4,2)]),
      4*p**2 + mmax([V('L',a,p) + V('S',b,p) + V('L',c,p) for a,b,c in compositions(n-4,3)]),
      3*p**2 + mmax([V('L',a,p) + V('L',b,p) for a,b in compositions(n-4,2)]),
      2*p**2 + V('L',n-4,p)
    ])


def R(n, p):
  return V('R', n, p) / OPT(n, p)
