
import numpy as np






def return_shuffled_inds(ind1,ind2):
    SHUFFLED = True

    while SHUFFLED:
        indices_rand1 = ind1
        indices_rand2 = ind2
        np.random.shuffle(indices_rand1)
        np.random.shuffle(indices_rand2)
        if len(np.where(indices_rand1-indices_rand2==0)[0])==0:
            SHUFFLED=False
    return ind1,ind2