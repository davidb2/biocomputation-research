#!/usr/bin/env python3.6
import argparse
import numpy as np

from decimal import Decimal
from functools import lru_cache

@lru_cache(maxsize=None)
def E(n)


def main(args):
  N = args.N
  loop = args.loop

  if loop:
    for i in range(1,N+1):
      a = compute(i)
      print(i, a, a/i)
  else:
    a = compute(N)
    b = (3*N-8-(-2)**(3-N))/18
    print(N, a, a/N)
    print(b)

def compute(N):
  M = np.zeros((N, N), dtype=np.float64)
  b = np.zeros((N, 1), dtype=np.float64)

  for i in range(N):
    b[i] = -sum(
      1/2**(j+1) for j in range(1,N+1) if i+2*j+1< N
    )

    M[i,i] = -1.
    for j in range(N):
      if i+2*j+1>=N: break
      M[i,i+2*j+1] = 1/2**(j+1)

  # print(M)
  print(b)
  # input()
  x = np.linalg.solve(M,b)
  return x[0,0]


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-N', type=int, required=True)
  parser.add_argument('-l', '--loop', action='store_true')
  args = parser.parse_args()
  main(args)
