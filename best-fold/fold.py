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

P = None
M = np.array([
  [.5, .5, 0.],
  [.5, .5, 0.],
  [.5, .5, 0.],    # <---- Initial state.
])

class Direction(Enum):
  RIGHT = +1
  LEFT = -1

@lru_cache(maxsize=None)
def p(i, j, delta):
  '''
  Computes Pr[X_delta = j | X_0 = i].
  '''
  assert delta >= 0
  if delta == 0: return float(i == j)
  return sum(M[i,k]*p(k,j,delta-1) for k in (0,1,2))

@lru_cache(maxsize=None)
def match(a, b):
  if a == b-1 or b == a-1: return 0
  #return int(P[a] == P[b] == 1)
  return p(2, 1, a)*p(1, 1, b-a)


@lru_cache(maxsize=None)
def score(i, j, k, kk):
  '''
  (Right)                              (Left)
  i -------> k                         k <------ i
             |           or            |
             v                         v
     j <--- kk                        kk ------------> j
  '''
  if kk > j or k < i: return 0
  return match(k, kk) + score(i, j, k-1, kk+1)


@lru_cache(maxsize=None)
def fold(i, k):
  N = len(P)
  if i == k or k == N-1: return 0
  return max(
    fold(k+1, j) + score(i, j, k, k+1)
    for j in range(k+1, N)
  )


def foldLoop(P):
  N = len(P)
  FOLD = np.zeros((N, N), dtype=np.float32)
  nextZag = np.zeros((N, N, 2), dtype=np.int32)

  # Explicit base case.
  FOLD[:, N-1] = 0
  for x in range(N):
    FOLD[x, x] = 0

  for i in reversed(range(N)):
    for k in reversed(range(i+1, N-1)):
      scores = [
        score(i, j, k, k+1) + FOLD[k+1, j] if j > k else 0
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


def plotZigZag(P, folds):
  edges = []
  bonds = []
  aacids = []

  fig, ax = plt.subplots(num=None, figsize=(32, 24), dpi=120, facecolor='w', edgecolor='k')
  def draw(a, b, i):
    color = ['blue', 'blue'][P[i]]
    ax.plot(*zip(a, b), color='black', linestyle='-', zorder=0)
    ax.plot(*zip(a, b), color=color, marker='o', linestyle='', markersize=10, zorder=1)

  N = len(P)
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
    if P[j] == 0: continue

    down = tuple(position + np.array([+1, 0]))
    right = tuple(position + np.array([0, +1]))

    linkedNeighbors = set()
    if j > 0:
      linkedNeighbors.add(tuple(positionsList[j-1]))
    if j < N-1:
      linkedNeighbors.add(tuple(positionsList[j+1]))

    for k, neighbor in enumerate([down, right]):
      if neighbor not in positionsSet or P[positionsSet[neighbor]] == 0: continue

      if neighbor not in linkedNeighbors and neighbor in positionsSet:
        xs, ys = zip(position, neighbor)

        # Should only happen for right turns.
        # if k == 0:
        #   xx = np.linspace(0, 1, 1000)
        #   yy = 0.1 * np.sin(32 * (2 * np.pi * xx))
        #   ax.plot(
        #     position[0] + xx,
        #     position[1] + yy,
        #     color='green', linestyle='-', zorder=0, alpha=0.25,
        #   )

  maxX, maxY = np.max(np.array(positionsList).T, axis=1)
  m = max(maxX, maxY)
  ax.set_xticks(np.linspace(0, m, m+1))
  ax.set_yticks(np.linspace(0, m, m+1))
  ax.axis('off')
  fig.canvas.set_window_title('Zig Zag Fold')
  plt.show()


def main(args):
  global P
  P = args.P

  ans = max(fold(0, k) for k in range(len(P)))
  print(ans)

  ans2, folds = foldLoop(P)
  print(ans2, folds)
  plotZigZag(P, folds)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-P', type=str, required=True)
  args = parser.parse_args()

  args.P = [int(x) for x in args.P]
  main(args)
