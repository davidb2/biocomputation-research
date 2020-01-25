#!/usr/bin/env python3.6
import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import sympy

from collections import namedtuple

# - `name` stores its symbol.
class NonTerminal:
  def __init__(self, name):
    assert len(name) == 1, len(name)
    self.name = name

  def __repr__(self):
    return self.name

# - `name` stores its symbol.
# - `values` is a list of possible values.
# - `pmf` maps a value to a probability (mass).
#   - Note that this scheme assumes independence.
class PTerminal:
  def __init__(self, name, values, pmf):
    assert len(name) == 1, len(name)
    self.name = name
    self.values = values
    self.pmf = pmf

  def __repr__(self):
    return self.name

# - `nonTerminal` is the source.
# - `rule` is a list of `Terminal`s and `PTerminal`s corresponding to
class Production:
  def __init__(self, nonTerminal, rule):
    self.nonTerminal = nonTerminal
    self.rule = rule

  # def score(self, values):
  #   pts = [pt for pt in self.rule if type(pt) is PTerminal]
  #   assert len(values) == len(pts), (len(values), len(pts))
  #   return np.prod([pt.pmf(value) for pt, value in zip(pts, values)])

  def __repr__(self):
    return (self.nonTerminal, self.rule)

# - `nonTerminal` is the source.
# - `productions` is a list `Production`s.
# Note: Undefined behavior if inputs are changed after instantiation.
class PCFG:
  def __init__(self, pTerminals, nonTerminals, productions=[]):
    names = set([pnt.name for pnt in pTerminals + nonTerminals])
    # Avoid name conflicts.
    assert len(names) == len(pTerminals) + len(nonTerminals), len(names)

    self.pTerminals = pTerminals
    self.nonTerminals = nonTerminals

    self.nameMap = {pnt.name: pnt for pnt in pTerminals + nonTerminals}

    self.productions = []
    for production in productions:
      self.addProduction(production)

  def addProduction(self, production):
    assert production.nonTerminal in self.nameMap, production.nonTerminal

    rule = []
    for c in production.rule:
      assert (c in self.nameMap), c
      rule.append(self.nameMap[c])

    self.productions.append(Production(
      nonTerminal=NonTerminal(production.nonTerminal),
      rule=rule,
    ))


pTerminals = [
  PTerminal(
    name='x',
    values=['0', '1'],
    pmf=lambda x: 0.5 if x in ['0', '1'] else 0.0,
  ),
]

nonTerminals = [
  NonTerminal('T'),
]

productions = [
  Production('T', 'xTx'),
]

pcfg = PCFG(pTerminals, nonTerminals, productions)
