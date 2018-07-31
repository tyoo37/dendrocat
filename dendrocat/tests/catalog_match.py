from astropy.table import Table
import dendrocat

table = Table.read('/users/bmcclell/dendrocat/dendrocat/data/EVLA_VLA_PointSourcePhotometry.ipac', format='ipac')

# load mc


mc.catalog_match(table, ra_colname='gracen', dec_colname='gdeccen')

