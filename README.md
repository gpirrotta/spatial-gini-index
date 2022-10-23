## Spatial Gini Index
A Python library to compute the **Spatial Gini Index**, a new indicator to measure spatial concentration.


Requirements
------------
* Python >=3.9


Installation
-------------

```pip install git+https://github.com/gpirrotta/spatial-gini-index.git```


Library Reference
-----------------

#### Class: sgi.metrics.SGI

**Constructor Parameters**
- **points**: `ndarray`, Nx2 array containing longitude/latitude data with `float` type
- **y**: `ndarray`, containing target values
- **step**: `str`, step distance method. Available `max-min` and `constant`. Default: `max-min`
- **label**: `ndarray`, points labels
- **permutation**: `int`, if different from 0 it use the permutation number to compute the p-value. Default: 0

**Methods**

* **sgi**: `float`, the Spatial Gini Index
* **target**: `ndarray`, the target values
* **pvalue_mc**: `float`, the pvalue
* **links**: `list`, the number of connections
* **variabilities**: `list`, the contiguity variability values
* **h_distances**: `list`, the h-distances values
* **distance_matrix**: 2D `ndarray` (NxN), the distance matrix
* **report**: pandas `dataframe`, general statistics
* **report_pvalue**: pandas `dataframe`, pvalue statistics


Using the tool
--------------

```python
import geopandas as gpd
import numpy as np
from sgi.metrics import SGI

my_territory_gdf = gpd.read_file('myterritory.shp')
centroids = my_territory_gdf.centroid
points = np.vstack([centroids.x, centroids.y]).T
labels = my_territory_gdf['TERRITORY_NAMES'].values

target = my_territory_gdf['MY_TARGET'].values

sgi = SGI(points, y, step=step, labels=labels)

print(sgi.sgi)
# 0.4739720213724825
            
print(sgi.report) # SEE TABLE BELOW
```
|   Spatial Lag |       h-distance |   Links |   Contiguity Variability | h-distance (links_name)   |     % Links |   % Contiguity Variability |   Cumulative Links |   Cumulative Contiguity Variability |   tg (rad) |   tg (degree) | SGI                | Time in seconds     |
|--------------:|-----------------:|--------:|-------------------------:|:--------------------------|------------:|---------------------------:|-------------------:|------------------------------------:|-----------:|--------------:|:-------------------|:--------------------|
|             1 |  77962.2         |     516 |              4.95042e+11 | Palermo-Trapani           | 0.0454946   |                0.0585515   |          0.0454946 |                           0.0585515 |   0.910237 |       52.1527 | 0.4739720213724825 | 0.08539298999999967 |
|             2 | 116740           |     554 |              4.1967e+11  | Palermo-Catania           | 0.048845    |                0.0496368   |          0.0943396 |                           0.108188  |   0.853671 |       48.9118 |                    |                     |
|             3 | 375746           |    4140 |              3.14259e+12 | Grosseto-Oristano         | 0.365015    |                0.371693    |          0.459355  |                           0.479881  |   0.807249 |       46.252  |                    |                     |
|             4 | 417920           |     596 |              4.83361e+11 | Palermo-Cagliari          | 0.0525481   |                0.05717     |          0.511903  |                           0.537051  |   0.809368 |       46.3734 |                    |                     |
|             5 | 454580           |     510 |              4.48709e+11 | Chieti-Messina            | 0.0449656   |                0.0530715   |          0.556868  |                           0.590122  |   0.814383 |       46.6607 |                    |                     |
|             6 | 492565           |     518 |              5.97971e+11 | Latina-Cremona            | 0.045671    |                0.0707256   |          0.602539  |                           0.660848  |   0.831518 |       47.6425 |                    |                     |
|             7 | 528304           |     452 |              4.97848e+11 | Pistoia-Sud Sardegna      | 0.0398519   |                0.0588834   |          0.642391  |                           0.719731  |   0.842117 |       48.2497 |                    |                     |
|             8 | 570481           |     528 |              3.18898e+11 | Crotone-Perugia           | 0.0465526   |                0.037718    |          0.688944  |                           0.757449  |   0.832726 |       47.7117 |                    |                     |
|             9 | 618608           |     522 |              3.38666e+11 | Salerno-La Spezia         | 0.0460236   |                0.040056    |          0.734967  |                           0.797505  |   0.826184 |       47.3369 |                    |                     |
|            10 | 666715           |     490 |              2.26715e+11 | Forli'-Cesena-Lecce       | 0.0432023   |                0.0268149   |          0.77817   |                           0.82432   |   0.81419  |       46.6496 |                    |                     |
|            11 | 712725           |     424 |              2.40208e+11 | Grosseto-Siracusa         | 0.0373832   |                0.0284108   |          0.815553  |                           0.852731  |   0.80768  |       46.2766 |                    |                     |
|            12 | 759239           |     342 |              1.44105e+11 | Salerno-Asti              | 0.0301534   |                0.0170441   |          0.845706  |                           0.869775  |   0.799428 |       45.8038 |                    |                     |
|            13 | 804866           |     302 |              1.32443e+11 | Pisa-Siracusa             | 0.0266267   |                0.0156648   |          0.872333  |                           0.88544   |   0.792855 |       45.4272 |                    |                     |
|            14 | 850903           |     292 |              9.67944e+10 | Ravenna-Siracusa          | 0.025745    |                0.0114484   |          0.898078  |                           0.896888  |   0.784735 |       44.962  |                    |                     |
|            15 | 899838           |     280 |              1.72958e+11 | Modena-Siracusa           | 0.024687    |                0.0204568   |          0.922765  |                           0.917345  |   0.782453 |       44.8312 |                    |                     |
|            16 | 942642           |     268 |              1.32279e+11 | Parma-Siracusa            | 0.023629    |                0.0156455   |          0.946394  |                           0.932991  |   0.778267 |       44.5914 |                    |                     |
|            17 | 985851           |     238 |              2.98013e+11 | Verona-Siracusa           | 0.020984    |                0.0352477   |          0.967378  |                           0.968238  |   0.785843 |       45.0255 |                    |                     |
|            18 |      1.03666e+06 |     222 |              1.67339e+11 | Brescia-Siracusa          | 0.0195733   |                0.0197922   |          0.986951  |                           0.988031  |   0.785945 |       45.0313 |                    |                     |
|            19 |      1.08797e+06 |     102 |              8.86412e+10 | Lecco-Siracusa            | 0.00899312  |                0.0104841   |          0.995944  |                           0.998515  |   0.786687 |       45.0738 |                    |                     |
|            20 |      1.11111e+06 |      34 |              1.23607e+10 | Bolzano-Siracusa          | 0.00299771  |                0.00146198  |          0.998942  |                           0.999977  |   0.785916 |       45.0297 |                    |                     |
|            21 |      1.15371e+06 |      10 |              1.78355e+08 | Aosta-Siracusa            | 0.000881679 |                2.10951e-05 |          0.999824  |                           0.999998  |   0.785485 |       45.005  |                    |                     |
|            22 |      1.15403e+06 |       2 |              1.87394e+07 | Aosta-Ragusa              | 0.000176336 |                2.21642e-06 |          1         |                           1         |   0.785398 |       45      |                    |



See the paper (Not ready yet)


License
-------
[MIT license](https://spdx.org/licenses/MIT.html)


Credits
-------
[Federico Benassi](mailto:benassi@istat.it)
[Massimo Mucciardi](mailto:massimo.mucciardi@unime.it)
[Giovanni Pirrotta](mailto:giovanni.pirrotta@unime.it)
