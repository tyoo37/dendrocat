import numpy as np

def mask(reg, cutout):
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')
    
def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5
