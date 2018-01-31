import numpy as np
import healpy
import aipy as a
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


flux = 0.0
frequency = 750 # MHz
catalog = 'misc'
# catalog = 'nvss'
# catalog = 'wenss'
src = '%f/%f' % (flux, frequency / 1.0e3)
srclist, cutoff, catalogs = a.scripting.parse_srcs(src, catalog)
cat = a.src.get_catalog(srclist, cutoff, catalogs)
nsrc = len(cat) # number of sources in cat
# print help(cat)
# print help(cat.values()[0])
# print [ str(cat.values()[i]) for i in range(nsrc) ]
# print str(cat.values()[0]).split()
# print [ cat.values()[i].name for i in range(nsrc) ]
valid_inds = []
for i in range(nsrc):
    if len(str(cat.values()[i]).split()) > 2:
        valid_inds.append(i)

# select sources
names = [ str(cat.values()[i]).split()[0] for i in valid_inds ]
ras = [ np.degrees(cat.values()[i]._ra) for i in valid_inds ]
decs = [ np.degrees(cat.values()[i]._dec) for i in valid_inds ]
jys = [ cat.values()[i].get_jys() for i in valid_inds ]

ras = np.array(ras)
decs = np.array(decs)
jys = np.array(jys)


# lon_0 is central longitude of projection.
m = Basemap(projection='moll', lon_0=0, celestial=True)

# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,30.))

# plt.title("Mollweide Projection")
# convert to map projection coords.
# Note that ra, dec can be scalars, lists or numpy arrays.
xpt, ypt = m(ras, decs)
# convert back to ra/dec
rapt, decpt = m(xpt, ypt, inverse=True)
# m.plot(xpt, ypt, 'bo')  # plot a blue dot there
# m.scatter(xpt, ypt, s=jys, color='r', alpha=.9)
m.scatter(xpt, ypt, s=jys/10, c='#9a0200', edgecolors='k', alpha=1) # deep red
# put some text next to the dot, offset a little bit
# (the offset is in map projection coordinates)
# plt.text(xpt+100000, ypt+100000, 'Boulder (%5.1fW,%3.1fN)' % (rapt, decpt))
for i, nm in enumerate(names):
    plt.text(xpt[i]-1000000, ypt[i]+500000, nm)
plt.savefig('10_sources.png')
plt.close()
