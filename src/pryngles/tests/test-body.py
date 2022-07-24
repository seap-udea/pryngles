##################################################################
#                                                                #
#.#####...#####...##..##..##..##...####...##......######...####..#
#.##..##..##..##...####...###.##..##......##......##......##.....#
#.#####...#####.....##....##.###..##.###..##......####.....####..#
#.##......##..##....##....##..##..##..##..##......##..........##.#
#.##......##..##....##....##..##...####...######..######...####..#
#................................................................#
#                                                                #
# PlanetaRY spanGLES                                             #
# The bright-side of the light-curve of (ringed) exoplanets      #
#                                                                #
##################################################################
# Jorge I. Zuluaga, Mario Sucerquia, Jaime A. Alvarado (C) 2022  #
##################################################################
import unittest
from pryngles import *
class Test(unittest.TestCase):
    def test_body(self):
        class defaults(object):
            orbit=dict()
            physics=dict()
            optics=dict()
        P=Body(defaults,"Test",None,dict(),dict(),dict())
        obj=Body(defaults,"Test",P,dict(),dict(),dict())
        obj._update_parent("parent")
        obj._update_childs("child1")
        obj._update_childs("child2")
        print(obj.childs)
        print(obj.hash)
        self.assertEqual([obj.parent],["parent"],True)
        self.assertEqual(obj.childs,["child1","child2"],True)
        self.assertRaises(AssertionError,
                          lambda:Body(defaults,"Test",1,dict(),dict(),dict())
                         )   


if __name__=="__main__":
        unittest.main(argv=['first-arg-is-ignored'],exit=False)
