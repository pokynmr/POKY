
import __main__
s = __main__.main_session
proj = s.project

#
# This is an example script for simple 2D metabolomics by boxing
# integration, PCA and clustering. Plots PCA and dendogram.
# ROIs (region of interests) are multivariate features.
#
# Developed by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
# Last update: January 17, 2021
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module
#
# Method:
#     Choose one reference spectrum to use,
#     Use "vd" to duplicate views and specify ROIs
#

########################################
# USER PARAMETER START
# Set a reference spectrum that has multiple views (made by vd)
# for ROIs. The view with the spectrum name will be excluded.
# The user should adjust view window region which is *ROI*
ref_specname = "reference spectrum name"

# Describe spectra to use.
specname_list = ['spectrum name 1',
                 'spectrum name 2',
                 'spectrum name 3',
                 'spectrum name 4',]
# If you wish to use all spectra, uncomment line below.
#specname_list = list(map(lambda sp: sp.name, proj.spectrum_list()))

# You can also exclude spectra from reference library using 
# excl_condname_list in case they are loaded.

# You can choose either, manual ROI list set up,
# or automatic set up using view windows (make_regions() function)

# Manual ROI list (uncomment and define to activate).
#ROI_list = [ ['metabolite 1', 'view name 1'],
#             ['metabolite 2', 'view name 2'],
#             ['metabolite 3', 'view name 3'],
#             ['metabolite 4', 'view name 4'],

# Condition for consideration (Can be set by two-letter-code "st"). 
# All condition will be considered if not set.
# condname_list = ['condition name']
# Condition to disregard. Such as BMRB_Standard_Library
# excl_condname_list = ['BMRB_Standard_Library',]

max_area = 10000 # maximum allowed area size for a ROI
# USER PARAMETER END
###########################################

from sputil import name_to_spectrum, name_to_view
# Establish condition list if not set
try:
  cond_list = list(map(lambda x: name_to_condition(x, s),
                    condname_list))
except:
  cond_list = list(map(lambda x: x, s.proj.condition_list())) 

try:
  excl_cond_list = list(map(lambda x: name_to_condition(x, s),
                    excl_condname_list))
except:
  excl_cond_list = []

# Exclude from cond_list using excl_cond_list
cond_list = [x for x in cond_list if x not in excl_cond_list]

# define function for making ROI list
# this can be replaced by custom definition like above
# we can use view name. we do not cluster by view names here
# however, it can be considered to make a rule to do that.
def make_ROIs(spec):
  ROI_list = []
  for view in s.project.view_list():
    if view.spectrum == spec and view.name != spec.name:
      ROI_list.append([view.name, view.name])
  return ROI_list

# for unit conversion
ref_spec = name_to_spectrum(ref_specname, s)
ppm_per_pt = (ref_spec.spectrum_width[0] / (ref_spec.data_size[0] - 1),
              ref_spec.spectrum_width[1] / (ref_spec.data_size[1] - 1))

# define function for boxing integration
# spec is a spectrum to analyze
from itertools import product
def box_sum(spec, region):
  npoints = (int((region[1][0] - region[0][0]) / ppm_per_pt[0]),
            int((region[1][1] - region[0][1]) / ppm_per_pt[1]))
  
  it_list = list(product(range(npoints[0]), range(npoints[1])))
  hts = []
  for i in range(npoints[0]):
    x = region[1][0] - ppm_per_pt[0] * i
    for j in range(npoints[1]):
      y = region[1][1] - ppm_per_pt[1] * j
      hts.append(spec.data_height((x, y)))
  sum_hts = sum(hts)
  avg_hts = sum_hts / (npoints[0] * npoints[1])
  return sum_hts, avg_hts

# Check if the user defined ROIs.
try:
  test = ROI_list[0]
except:
  ROI_list = make_ROIs(ref_spec)

spec_list = list(map(lambda sp: name_to_spectrum(sp, s), 
                specname_list))

# Filter spectrum by condition
filt_list =list(map(lambda x: x.condition in cond_list, spec_list))

view_list = list(map(lambda v: name_to_view(v[1], s), ROI_list))
name_list = list(map(lambda v: v[0], ROI_list))

nspec = len(spec_list) # Num of examples
nfeat = len(view_list) # Num of features

# Check if any view area > max_area
area_list = list(map(lambda v: 
        (v.region[1][0] - v.region[0][0]) * \
        (v.region[1][1] - v.region[0][1]), 
        view_list))
sorted_area_list = sorted(area_list)
if sorted_area_list[0] > max_area:
  print('Maximum area error. Make sure your ROIs correctly defined.')
  raise SystemExit

# Import data science modules
from sklearn.preprocessing import StandardScaler 
from sklearn.decomposition import PCA
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

# Calculate ROI box integration
# Make a 2D-matrix array and save in dataframe
hts_list = []
for spec in spec_list:
  hts_list.append([])
  for i in range(len(view_list)):
    v = view_list[i]
    sum_hts, avg_hts = box_sum(spec, v.region)
    hts_list[-1].append(sum_hts)

df = pd.DataFrame(hts_list, columns = name_list, index=specname_list, 
                            dtype=float)

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

print('\n* Box sum values')
print(df.head())

# Normalization. Should we need to apply center averaging before this?
norm_features = StandardScaler().fit_transform(df.values)
norm_features_df = pd.DataFrame(norm_features, columns=df.columns,
                                index=df.index)
print('\n* Normalized box sum values')
print(norm_features_df.head())

pca = PCA(n_components=2)
pca_data = pca.fit_transform(norm_features)
print('\n* PCA results')
print(pca_data)

print('\n* Loadings')
loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], 
                        index=name_list)
print(loadings.head())

# Clustering
from sklearn.cluster import AgglomerativeClustering
cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean', 
                                  linkage='ward')
cluster.fit_predict(pca_data)
print('\n* Agglomerative Clustering')
print(cluster.labels_)

# Plotting
# Scatter PCA plot
plt.figure()
plt.scatter(pca_data[:,0], pca_data[:,1], c=cluster.labels_, 
            cmap='rainbow')
plt.scatter(loadings.value[:,0], loadings.value[:,1], c='gray')
for label, x, y in zip(loadings.index, loadings.value(:,0),
                       loadings.value(:,1)):
plt.annotate(label, 
             xy=(x, y), xytext=(3,3),
             textcoords='offset points', ha='right', va='bottom')
plt.pause(0.1)

# Dendogram
from scipy.cluster.hierarchy import dendrogram
children = cluster.children_
distance = np.arange(children.shape[0])
nobs = np.arange(2, children.shape[0]+2)
lmat = np.column_stack([children, distsance, nobs]).astyle(float)
dendogran(lmat, labels=cluster.labels_)
plt.figure()

plt.show(block=False)