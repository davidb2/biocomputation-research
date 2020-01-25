#!/usr/bin/env python3.6
import argparse
import numpy as np
import sys


def findFirst(A, i):
  for j in range(i+3, len(A), 2):
    if A[j] == 'H': return j
  return None

def greedy(A, i=0):
  if i == len(A): return 0
  if A[i] == 'P': return greedy(A, i+1)
  js = findFirst(A, i)
  if js is None: return 0
  return 1 + greedy(A, js)


def best(A):
  return 2*min((A[::2]=='H').sum(), (A[1::2]=='H').sum()) + 2

def ratio(A):
  return greedy(A)/best(A)

def main(args):
  N = args.N
  S = args.S

  sys.setrecursionlimit(N+10)
  A = np.random.choice(['H','P'], (S, N))
  # print(ratio(A[0]))
  print(A)

  ans = np.mean([ratio(X) for X in A])
  print(ans)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-N', type=int, required=True)
  parser.add_argument('-S', type=int, default=1000)

  args = parser.parse_args()
  main(args)
