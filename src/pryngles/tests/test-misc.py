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
#                                                                #
##################################################################
# License http://github.com/seap-udea/pryngles-public            #
##################################################################
import unittest
from pryngles import *
class Test(unittest.TestCase):
	def test_misc(self):
	
	    #Get path
	    filepath=Misc.get_data("diffuse_reflection_function.data")
	    print(filepath)
	
	    #print_df dataframe
	    import pandas as pd
	    import numpy as np
	    df=pd.DataFrame(np.zeros((5,3)),columns=["a","b","c"])
	    Misc.print_df(df)
	
	    #Flatten
	    print(list(Misc.flatten(["hola"])))
	    print(list(Misc.flatten(["hola",["perro","gato"]])))
	
	    #Get methods
	    print(Misc.get_methods(Misc))
	    
	    #Hash
	    d=dict(a=1,b=3,c=np)
	    print(Misc.calc_hash(d))
	    P=PrynglesCommon()
	    print(Misc.calc_hash(P))
	    print(Misc.calc_hash(PrynglesCommon))
	
	    Verbose.save_test_fig('misc-test_misc')
	    plt.close('all')

	def test_time(self):
	    Verbose.VERBOSITY=VERB_ALL
	    
	    print(Misc.TIME_IMPORT)
	    Misc.elapsed_time()
	
	    Misc.elapsed_time(show=False)
	    time.sleep(2)
	    Misc.elapsed_time(total=True)
	
	    Misc.elapsed_time(imp=True)
	    Verbose.VERBOSITY=VERB_NONE    
	    
	    Verbose.save_test_fig('misc-test_time')
	    plt.close('all')

if __name__=="__main__":
   unittest.main(argv=['first-arg-is-ignored'],exit=False)
