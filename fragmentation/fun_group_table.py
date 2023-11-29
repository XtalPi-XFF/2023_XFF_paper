#!/usr/bin/env python
# coding=utf-8
'''define functional groups and their scores'''
from __future__ import absolute_import
from base_fragment import BaseFragment, BaseAtom, BaseBond

#degree 1
basefrag0 = BaseFragment([BaseAtom('C', 0, score=0)], [])
basefrag1 = BaseFragment([BaseAtom('Y', 0, score=1)], [])

#degree 2
basefrag2_1 = BaseFragment([BaseAtom('C', 0, score=0), BaseAtom('C', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2)])
basefrag2_2 = BaseFragment([BaseAtom('C', 0, score=0), BaseAtom('Y', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2)])
basefrag3_1 = BaseFragment([BaseAtom('C', 0, score=7), BaseAtom('C', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=2)])
basefrag3_2 = BaseFragment([BaseAtom('C', 0, score=7), BaseAtom('Y', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=2)])
basefrag3_3 = BaseFragment([BaseAtom('C', 0, score=7), BaseAtom('C', 1), BaseAtom('Y', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=2)])
basefrag3_4 = BaseFragment([BaseAtom('C', 0, score=7), BaseAtom('Y', 1), BaseAtom('Y', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=2)])
basefrag3_5 = BaseFragment([BaseAtom('Y', 0, score=7), BaseAtom('C', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=2)])
basefrag3_6 = BaseFragment([BaseAtom('Y', 0, score=7), BaseAtom('C', 1), BaseAtom('Y', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=2)])
basefrag3_7 = BaseFragment([BaseAtom('Y', 0, score=7), BaseAtom('Y', 1), BaseAtom('Y', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=2)])
basefrag4_1 = BaseFragment([BaseAtom('C', 0, score=10), BaseAtom('C', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=3)])
basefrag4_2 = BaseFragment([BaseAtom('C', 0, score=10), BaseAtom('C', 1), BaseAtom('Y', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=3)])
basefrag4_3 = BaseFragment([BaseAtom('C', 0, score=10), BaseAtom('Y', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=3)])
basefrag4_4 = BaseFragment([BaseAtom('C', 0, score=10), BaseAtom('Y', 1), BaseAtom('Y', 2)], [BaseBond(0, 1), BaseBond(0, 2, bondtype=3)])
basefrag5 = BaseFragment([BaseAtom('Y', 0, score=1), BaseAtom('C', 1), BaseAtom('C', 2)], [BaseBond(0, 1), BaseBond(0, 2)])

#degree3
basefrag6_1 = BaseFragment([BaseAtom('C', 0, score=2), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(3, 4)])
basefrag6_2 = BaseFragment([BaseAtom('C', 0, score=4), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])
basefrag7_1 = BaseFragment([BaseAtom('C', 0, score=9), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1, bondtype=2), BaseBond(0, 2), BaseBond(0, 3), BaseBond(3, 4)])
basefrag7_2 = BaseFragment([BaseAtom('C', 0, score=8), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3)],
                           [BaseBond(0, 1, bondtype=2), BaseBond(0, 2), BaseBond(0, 3)])
basefrag8_1 = BaseFragment([BaseAtom('C', 0, score=11), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('Y', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(3, 4)])
basefrag8_2 = BaseFragment([BaseAtom('C', 0, score=15), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])
basefrag8_3 = BaseFragment([BaseAtom('C', 0, score=15), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])
basefrag8_4 = BaseFragment([BaseAtom('C', 0, score=15), BaseAtom('Y', 1), BaseAtom('Y', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])

basefrag9_1 = BaseFragment([BaseAtom('Y', 0, score=20), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(3, 4)])
basefrag9_2 = BaseFragment([BaseAtom('Y', 0, score=20), BaseAtom('Y', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(3, 4)])
basefrag9_3 = BaseFragment([BaseAtom('Y', 0, score=20), BaseAtom('Y', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(3, 4)])
basefrag9_4 = BaseFragment([BaseAtom('Y', 0, score=22), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])
basefrag10_1 = BaseFragment([BaseAtom('Y', 0, score=27), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('Y', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(3, 4)])
basefrag10_2 = BaseFragment([BaseAtom('Y', 0, score=31), BaseAtom('Y', 1), BaseAtom('Y', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])
basefrag10_3 = BaseFragment([BaseAtom('Y', 0, score=31), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])
basefrag10_4 = BaseFragment([BaseAtom('Y', 0, score=31), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3)])
basefrag11_1 = BaseFragment([BaseAtom('C', 0, score=26), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('Y', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3, bondtype=2), BaseBond(3, 4)])
basefrag11_2 = BaseFragment([BaseAtom('C', 0, score=25), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3, bondtype=2)])
basefrag11_3 = BaseFragment([BaseAtom('C', 0, score=25), BaseAtom('Y', 1), BaseAtom('C', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3, bondtype=2)])
basefrag11_4 = BaseFragment([BaseAtom('C', 0, score=25), BaseAtom('Y', 1), BaseAtom('Y', 2), BaseAtom('Y', 3)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3, bondtype=2)])

#degree4
basefrag12_1 = BaseFragment([BaseAtom('C', 0, score=3), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4), BaseAtom('C', 5), BaseAtom('C', 6)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5), BaseBond(2, 6)])
basefrag12_2 = BaseFragment([BaseAtom('C', 0, score=5), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4), BaseAtom('C', 5)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(2, 5)])
basefrag12_3 = BaseFragment([BaseAtom('C', 0, score=6), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])
basefrag13_1 = BaseFragment([BaseAtom('C', 0, score=12), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5), BaseAtom('C', 6)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5), BaseBond(2, 6)])
basefrag13_2 = BaseFragment([BaseAtom('C', 0, score=14), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5)])
basefrag13_3 = BaseFragment([BaseAtom('C', 0, score=16), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(2, 5)])
basefrag13_4 = BaseFragment([BaseAtom('C', 0, score=18), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4)],
                           [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])
basefrag14_1 = BaseFragment([BaseAtom('C', 0, score=13), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5), BaseAtom('C', 6)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5), BaseBond(2, 6)])
basefrag14_2 = BaseFragment([BaseAtom('C', 0, score=17), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(2, 5)])
basefrag14_3 = BaseFragment([BaseAtom('C', 0, score=19), BaseAtom('Y', 1), BaseAtom('Y', 2), BaseAtom('Y', 3), BaseAtom('Y', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])
basefrag14_4 = BaseFragment([BaseAtom('C', 0, score=19), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('Y', 3), BaseAtom('Y', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])
basefrag14_5 = BaseFragment([BaseAtom('C', 0, score=19), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('Y', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])

basefrag15_1 = BaseFragment([BaseAtom('Y', 0, score=21), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4), BaseAtom('C', 5), BaseAtom('C', 6)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5), BaseBond(2, 6)])
basefrag15_2 = BaseFragment([BaseAtom('Y', 0, score=23), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4), BaseAtom('C', 5)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(2, 5)])
basefrag15_3 = BaseFragment([BaseAtom('Y', 0, score=24), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('C', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])
basefrag16_1 = BaseFragment([BaseAtom('Y', 0, score=28), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5), BaseAtom('C', 6)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5), BaseBond(2, 6)])
basefrag16_2 = BaseFragment([BaseAtom('Y', 0, score=30), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5)])
basefrag16_3 = BaseFragment([BaseAtom('Y', 0, score=32), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(2, 5)])
basefrag16_4 = BaseFragment([BaseAtom('Y', 0, score=34), BaseAtom('C', 1), BaseAtom('C', 2), BaseAtom('C', 3), BaseAtom('Y', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])

basefrag17_1 = BaseFragment([BaseAtom('Y', 0, score=29), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5), BaseAtom('C', 6)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(4, 5), BaseBond(2, 6)])
basefrag17_2 = BaseFragment([BaseAtom('Y', 0, score=33), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('Y', 4), BaseAtom('C', 5)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4), BaseBond(2, 5)])
basefrag17_3 = BaseFragment([BaseAtom('Y', 0, score=35), BaseAtom('C', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('Y', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])
basefrag17_4 = BaseFragment([BaseAtom('Y', 0, score=35), BaseAtom('Y', 1), BaseAtom('Y', 2), BaseAtom('C', 3), BaseAtom('Y', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])
basefrag17_5 = BaseFragment([BaseAtom('Y', 0, score=35), BaseAtom('Y', 1), BaseAtom('Y', 2), BaseAtom('Y', 3), BaseAtom('Y', 4)],
                            [BaseBond(0, 1), BaseBond(0, 2), BaseBond(0, 3), BaseBond(0, 4)])

#degree gt4


fun_grop_dict = {
    1:[basefrag0, basefrag1],
    2:[basefrag3_1, basefrag3_2, basefrag3_3, basefrag3_4, basefrag3_5, basefrag3_6, basefrag3_7, basefrag4_1, basefrag4_2, basefrag4_3, basefrag4_4, basefrag5, basefrag2_2, basefrag2_1],
    3:[basefrag6_1, basefrag6_2, basefrag7_1, basefrag7_2, basefrag9_1, basefrag9_2, basefrag9_3, basefrag9_4, basefrag10_1, basefrag10_2, basefrag10_3, basefrag10_4, basefrag11_1, basefrag11_2, basefrag11_3, basefrag11_4, basefrag8_1, basefrag8_2,basefrag8_3, basefrag8_4],
    4:[basefrag12_1, basefrag12_2, basefrag12_3, basefrag13_1, basefrag13_2, basefrag13_3, basefrag13_4, basefrag14_1, basefrag14_2, basefrag14_3, basefrag14_4, basefrag14_5, basefrag15_1, basefrag15_2, basefrag15_3, basefrag16_1, basefrag16_2, basefrag16_3, basefrag16_4, basefrag17_1, basefrag17_2, basefrag17_3, basefrag17_4, basefrag17_5]
}

#bigring
bigring_bond_score = {
    ('C', 'C', 2, 2):0,
    ('C', 'C', 3, 2):1,
    ('C', 'C', 2, 3):1,
    ('C', 'C', 2, 4):1,
    ('C', 'C', 4, 2):1,
    ('C', 'C', 3, 4):2,
    ('C', 'C', 4, 3):2,
    ('C', 'C', 4, 4):5,
    ('C', 'Y', 2, 2):2,
    ('Y', 'C', 2, 2):2,
    ('C', 'C', 3, 3):3,
    ('C', 'Y', 3, 2):4,
    ('C', 'Y', 2, 3):4,
    ('Y', 'C', 3, 2):4,
    ('Y', 'C', 2, 3):4,
    ('Y', 'C', 2, 4):4,
    ('Y', 'C', 4, 2):4,
    ('C', 'Y', 4, 2):4,
    ('C', 'Y', 2, 4):4,
    ('C', 'Y', 3, 3):5,
    ('Y', 'C', 3, 3):5,
    ('C', 'Y', 3, 4):6,
    ('Y', 'C', 3, 4):6,
    ('C', 'Y', 4, 3):6,
    ('Y', 'C', 4, 3):6,
    ('Y', 'Y', 2, 2):6,
    ('C', 'Y', 4, 4):7,
    ('Y', 'C', 4, 4):7,
    ('Y', 'Y', 2, 3):7,
    ('Y', 'Y', 2, 4):7,
    ('Y', 'Y', 4, 2):7,
    ('Y', 'Y', 3, 2):7,
    ('Y', 'Y', 3, 3):8,
    ('Y', 'Y', 3, 4):8,
    ('Y', 'Y', 4, 3):8,
    ('Y', 'Y', 4, 4):9
}
