


class MasterCatalog:
    """
    An object to store combined data from two or more RadioSource objects.
    """
    
    def __init__(self, *args, catalog=None):
        """
        Create a new master catalog object.
        
        Parameters
        ----------
        catalog : astropy.table.Table object
            The master table from which to build the catalog object.
        *args : radiosource.RadioSource objects
            RadioSource objects from which the master table was built.
        
        NOTE: Currently only supports multiple objects if they're of different
              frequencies (the freq_id must be unique).
        """
        if catalog is not None:
            self.catalog = catalog
        
        obj_prefix = 'radiosource_'
        
        for obj in args:
        
            if isinstance(obj, MasterCatalog):
                for key in obj.__dict__.keys():
                    if key.split('_')[0]+'_' == obj_prefix:
                        self.__dict__[key] = obj[key]
            else:        
                self.__dict__[obj_prefix+obj.freq_id] = obj
