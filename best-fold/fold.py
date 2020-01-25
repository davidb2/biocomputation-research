#!/usr/bin/env python3.6
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from enum import Enum
from functools import lru_cache
from matplotlib.textpath import TextPath

'''
Best accordian fold on a 2d square lattice.
'''
class Direction(Enum):
  RIGHT = +1
  LEFT = -1

class Fold:
  def __init__(self, mode, P=None, N=None, M=None):
    self.mode = mode
    self.N = N
    self.P = P
    self.M = M

    # Cache recursive functions per instance.
    self.p = lru_cache(maxsize=None)(self.p)
    self.match = lru_cache(maxsize=None)(self.match)
    self.score = lru_cache(maxsize=None)(self.score)
    self.fold = lru_cache(maxsize=None)(self.fold)

  def p(self, i, j, delta):
    '''
    Computes Pr[X_delta = j | X_0 = i].
    '''
    assert delta >= 0
    if delta == 0: return float(i == j)
    return sum(self.M[i,k]*self.p(k,j,delta-1) for k in (0,1,2))

  def match(self, a, b):
    if a == b-1 or b == a-1: return 0

    if self.mode == 'exact': return int(self.P[a] == self.P[b] == 1)
    if self.mode == 'random': return self.p(2, 1, a)*self.p(1, 1, b-a)


  def score(self, i, j, k, kk):
    '''
    (Right)                              (Left)
    i -------> k                         k <------ i
               |           or            |
               v                         v
       j <--- kk                        kk ------------> j
    '''
    if kk > j or k < i: return 0
    return self.match(k, kk) + self.score(i, j, k-1, kk+1)


  def fold(self, i, k):
    N = len(self)
    if i == k or k == N-1: return 0
    return max(
      self.fold(k+1, j) + self.score(i, j, k, k+1)
      for j in range(k+1, N)
    )

  def __len__(self):
    return len(self.P) if self.mode == 'exact' else self.N

  def foldLoop(self):
    N = len(self)
    FOLD = np.zeros((N, N), dtype=np.float32)
    nextZag = np.zeros((N, N, 2), dtype=np.int32)

    # Explicit base case.
    FOLD[:, N-1] = 0
    for x in range(N):
      FOLD[x, x] = 0

    for i in reversed(range(N)):
      for k in reversed(range(i+1, N-1)):
        scores = [
          self.score(i, j, k, k+1) + FOLD[k+1, j] if j > k else 0
          for j in range(N)
        ]

        # Find index of right-most maximum.
        jp = N - np.argmax(scores[::-1]) - 1
        FOLD[i, k] = scores[jp]
        nextZag[i, k] = (k+1, jp)

    # Find the best place to start folding.
    j = N - np.argmax(FOLD[0, ::-1]) - 1

    # Recontruct best fold.
    folds = [(0, j)]
    while folds[-1][-1] < N - 1:
      lastFold = folds[-1]
      folds.append(tuple(nextZag[lastFold[0], lastFold[1]]))

    return FOLD[0, j], folds


  def plotZigZag(self, folds):
    edges = []
    bonds = []
    aacids = []

    fig, ax = plt.subplots(num=None, figsize=(32, 24), dpi=120, facecolor='w', edgecolor='k')
    def draw(a, b, i):
      color = 'blue' if self.mode == 'random' or self.P[i] == 0 else 'orange'
      ax.plot(*zip(a, b), color='black', linestyle='-', zorder=0)
      ax.plot(*zip(a, b), color=color, marker='o', linestyle='', markersize=10, zorder=1)

    N = len(self)
    positionsSet = {}
    positionsList = []

    i = 0
    pos = np.array([0, 0])
    dc = np.array([0, +1])

    for k, fold in enumerate(folds):
      for _ in range(fold[1] - fold[0]):
        positionsSet[tuple(pos)] = i
        positionsList.append(pos)
        nextPos = pos + dc
        draw(pos, nextPos, i)
        pos = nextPos
        i += 1

      # Go down a row.
      positionsSet[tuple(pos)] = i
      positionsList.append(pos)

      if k == len(folds)-1: break

      nextPos = pos + np.array([+1, 0])
      draw(pos, nextPos, i)
      pos = nextPos
      dc *= -1
      i += 1

    # Draw bonds/connections.
    for j, position in enumerate(positionsList):
      if self.mode == 'exact' and self.P[j] == 0: continue

      down = tuple(position + np.array([+1, 0]))
      right = tuple(position + np.array([0, +1]))

      linkedNeighbors = set()
      if j > 0:
        linkedNeighbors.add(tuple(positionsList[j-1]))
      if j < N-1:
        linkedNeighbors.add(tuple(positionsList[j+1]))

      for k, neighbor in enumerate([down, right]):
        if neighbor not in positionsSet or (self.mode == 'exact' and self.P[positionsSet[neighbor]] == 0): continue

        if neighbor not in linkedNeighbors and neighbor in positionsSet:
          xs, ys = zip(position, neighbor)

          # Should only happen for right turns.
          if self.mode == 'exact' and k == 0:
            xx = np.linspace(0, 1, 1000)
            yy = 0.1 * np.sin(32 * (2 * np.pi * xx))
            ax.plot(
              position[0] + xx,
              position[1] + yy,
              color='green', linestyle='-', zorder=0, alpha=0.25,
            )

    maxX, maxY = np.max(np.array(positionsList).T, axis=1)
    m = max(maxX, maxY)
    ax.set_xticks(np.linspace(0, m, m+1))
    ax.set_yticks(np.linspace(0, m, m+1))
    ax.axis('off')
    fig.canvas.set_window_title('Zig Zag Fold')
    plt.show()


def main(args):
  F = Fold(mode=args.mode, N=args.N, P=args.P, M=np.array([
    [.5, .5, 0.],
    [.5, .5, 0.],
    [.5, .5, 0.],    # <---- Initial state.
  ]))

  ans = max(F.fold(0, k) for k in range(len(F)))
  print(ans)

  ans2, folds = F.foldLoop()
  print(ans2, folds)
  F.plotZigZag(folds)


if __name__ == '__main__':
  parser = argparse.ArgumentParser('Fold a polypeptide')
  subparsers = parser.add_subparsers(help='help for subcommands', dest='mode')
  subparsers.required = True

  parserExact = subparsers.add_parser('exact', help='Fold a given sequence of 0s and 1s')
  parserExact.add_argument('-P', type=str, required=True)

  parserRandom = subparsers.add_parser('random', help='Fold a random polypeptide of length N')
  parserRandom.add_argument('-N', type=int, required=True)

  args = parser.parse_args()

  args.N = args.__dict__.get('N')
  args.P = args.__dict__.get('P')
  args.P = [int(x) for x in args.P] if args.mode == 'exact' else args.P

  main(args)
