

# Keywords between dollar signs are substituted by subversion.
# The date and build will only tell the truth after a branch or tag,
# since different files in trunk will have been changed at different times
date ="$Date$".split()[1]
revision = "$Revision$".split()[1]


__version__ = '3.4' 


description = "CoreBio %s (%s)" % (__version__,  date)


